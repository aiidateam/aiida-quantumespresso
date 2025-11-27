"""Tools for nodes created by running the `PwCalculation` class."""
import numpy as np

from aiida.common import AttributeDict, exceptions
from aiida.tools.calculations.base import CalculationTools
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

from aiida_quantumespresso.calculations.functions.create_magnetic_configuration import create_magnetic_configuration

from xmlschema import XMLSchema
from xml.etree import ElementTree
from aiida_quantumespresso.parsers.parse_xml.parse import get_schema_filepath


class PwCalculationTools(CalculationTools):
    """Calculation tools for `PwCalculation`.

    Methods implemented here are available on any `CalcJobNode` produced by the `PwCalculation class through the `tools`
    attribute.
    """

    def get_scf_accuracy(self, index=0):
        """Return the array of SCF accuracy values for a given SCF cycle.

        :param index: the zero-based index of the desired SCF cycle
        :return: a list of SCF accuracy values of a certain SCF cycle.
        :raises ValueError: if the node does not have the `output_trajectory` output
        :raises ValueError: if `output_trajectory` does not have the `scf_accuracy` or `scf_iterations` arrays
        :raises IndexError: if the `index` is out of range
        """
        try:
            trajectory = self._node.outputs.output_trajectory
        except exceptions.NotExistent as exc:
            raise ValueError(f'{self._node} does not have the `output_trajectory` output node') from exc

        try:
            scf_accuracy = trajectory.get_array('scf_accuracy')
        except KeyError as exc:
            raise ValueError(f'{trajectory} does not contain the required `scf_accuracy` array') from exc

        try:
            scf_iterations = trajectory.get_array('scf_iterations')
        except KeyError as exc:
            raise ValueError(f'{trajectory} does not contain the required `scf_iterations` array') from exc

        number_of_frames = len(scf_iterations)

        if index < -number_of_frames or index >= number_of_frames:
            raise IndexError(f'invalid index {index}, must be between {-number_of_frames} and {number_of_frames - 1}')

        # building an scf_accuracy_index for easier manipulation
        scf_accuracy_index = [0]
        for i in scf_iterations:
            scf_accuracy_index.append(scf_accuracy_index[-1] + i)

        if index < 0:
            return scf_accuracy[scf_accuracy_index[index - 1] : scf_accuracy_index[index]]

        return scf_accuracy[scf_accuracy_index[index] : scf_accuracy_index[index + 1]]

    def get_magnetic_configuration(self, atol: float = 0.5, ztol: float = 0.05) -> AttributeDict:
        """Get the final magnetic configuration of a ``pw.x`` calculation."""
        try:
            structure = self._node.outputs.output_structure
        except AttributeError:
            structure = self._node.inputs.structure

        try:
            magnetic_moments = self._node.outputs.output_trajectory.get_array('atomic_magnetic_moments')[-1].tolist()
        except KeyError as exc:
            raise KeyError('Could not find the `atomic_magnetic_moments` in the `output_trajectory`.') from exc

        results = create_magnetic_configuration(
            structure,
            magnetic_moments,
            atol=atol,
            ztol=ztol,
            metadata={'store_provenance': False},
        )
        structure_kindname_position = [(site.kind_name, site.position) for site in structure.sites]
        allo_kindname_position = [(site.kind_name, site.position) for site in results['structure'].sites]
        structure_kindnames_sorted = sorted(structure_kindname_position, key=lambda el: el[1])
        allo_kindnames_sorted = sorted(allo_kindname_position, key=lambda el: el[1])

        requires_new_kinds = structure_kindnames_sorted != allo_kindnames_sorted
        non_magnetic = all(abs(magn) < ztol for magn in results['magnetic_moments'].get_dict().values())

        if requires_new_kinds:
            results = create_magnetic_configuration(
                structure,
                magnetic_moments,
                atol=atol,
                ztol=ztol,
                metadata={'store_provenance': False},
            )
            structure = results['results']

        return AttributeDict(
            {'structure': structure, 'magnetic_moments': None if non_magnetic else results['magnetic_moments']}
        )

    def get_occupations(self, reshape=False) -> list[dict]:
        """Return the occupations for a PwCalculation/PwBaseWorkChain node as a list of python dictionaries.
        Output depends on the type of calculation (unpolarized, collinear spin-polarized, non-collinear spin-polarized).

        Example of an entry in the returned list:
        
        Unpolarized:
        {
            'atom_index': 1,
            'kind_name': 'Fe',
            'manifold': 'd',
            'occupations': {
                'up-down': np.array([[...], [...], ...])
            }
        }
        Collinear spin-polarized:
        {
            'atom_index': 1,
            'kind_name': 'Fe',
            'manifold': 'd',
            'occupations': {
                'up': np.array([[...], [...], ...]),
                'down': np.array([[...], [...], ...])
            }
        }
        Non-collinear spin-polarized:
        {
            'atom_index': 1,
            'kind_name': 'Fe',
            'manifold': 'd',
            'occupations': {
                'up-down': np.array([[...], [...], ...])
            }
        }   
        
        """

        # assert first that this is a Hubbard calculation
        # so for instance if you try to call this method on a normal pw calculation it raises a ValueError
        if not isinstance(self._node.inputs.structure, HubbardStructureData):
             raise ValueError('The input structure is not a HubbardStructureData node.')

        
        with self._node.outputs.retrieved.base.repository.open('data-file-schema.xml') as xml_file:
            xml_parsed = ElementTree.parse(xml_file)

        schema_filepath = get_schema_filepath(xml_parsed)
        xsd = XMLSchema(schema_filepath)
        xml_dictionary = xsd.to_dict(xml_parsed, validation='skip')

        dft_u = xml_dictionary['output']['dft']['dftU']
        
        # We aggregate data by atom index
        aggregated_data = {}

        # Helper to standardize processing. 
        # We treat the two lists differently based on your specifications.
        # 1. Hubbard_ns_nc: Non-collinear (Single matrix, ignore @spin)
        # 2. Hubbard_ns:    Collinear (Split by @spin) OR No-Polarization (Single matrix, no @spin)
        
        source_lists = [
            (dft_u.get('Hubbard_ns_nc', []), True),  # (List, is_non_collinear)
            (dft_u.get('Hubbard_ns', []),    False)  # (List, is_standard_ns)
        ]

        # this needs to be verified together with the QE developers, but it seems that in the case of non-collinear
        # calculations, what is reported as "occupation_matrix" is actually the absolute value,
        # since it should have both real and imaginary parts

        for entries, is_non_collinear in source_lists:
            for entry in entries:
                atom_index = entry['@index']
                
                # The quantum espresso matrices are always as big as the biggest shell 
                # so if we have a d and a p orbital, the p shell matrix will also be 5x5 with zeros padding
                # same, if we have f and d, the d will be 7x7 with zeros padding
                # the @dims attribute for now does not reflect this
                # so we can for now extract a dimension based on the name of the manifold
                # 
                #
                # this will need to be done more intelligently (if we have Wannier orbitals there is no
                # guarantee for the dimension of the underlying shell) or QE should made to print the
                # correct dimensions in the @dims attribute
                actual_dim_map = {'s': 1, 'p': 3, 'd': 5, 'f': 7}


                # Parse matrix
                occ_matrix = np.array(entry['$'])

                # Initialize container if new
                if atom_index not in aggregated_data:
                    aggregated_data[atom_index] = {
                        'atom_index': atom_index,
                        'kind_name': entry['@specie'],
                        'manifold': entry['@label'],
                        'occupations': {} 
                    }
                # strip the number in the manifold label to get the actual dimension
                manifold_type = ''.join(filter(str.isalpha, entry['@label']))
                if manifold_type in actual_dim_map.keys():
                    actual_dim = actual_dim_map[manifold_type]
                    if is_non_collinear:

                        actual_dim *= 2  # Non-collinear has double size
                        # print(f"Non-collinear manifold detected, doubling actual_dim to {actual_dim}")
                else:
                    actual_dim = entry['@dims'][0]  # Fallback to reported dimension
                
                # reshape to @dims get the actual_dim x actual_dim block and reshape back to linear
                occ_matrix = occ_matrix.reshape(entry['@dims'])[:actual_dim, :actual_dim].reshape(-1)

                if reshape:
                    occ_matrix = occ_matrix.reshape((actual_dim, actual_dim))

                # --- LOGIC FOR THE 3 CASES ---

                # Case 3: Hubbard_ns_nc 
                # (@spin is always 1, but we treat it as a single full block)
                if is_non_collinear:
                    aggregated_data[atom_index]['occupations']['up-down'] = occ_matrix

                # Case 1 & 2: Hubbard_ns
                else:
                    spin_val = entry.get('@spin')
                    
                    if spin_val is None:
                        # Case 1: Hubbard_ns with no polarization (No @spin)
                        # Treat as single matrix
                        aggregated_data[atom_index]['occupations']['up-down'] = occ_matrix
                    else:
                        # Case 2: Hubbard_ns with polarization (@spin is 1 or 2)
                        # Treat as dictionary with 'up'/'down'
                        spin_label = 'up' if spin_val == 1 else 'down'
                             
                        aggregated_data[atom_index]['occupations'][spin_label] = occ_matrix

        return list(aggregated_data.values()) 