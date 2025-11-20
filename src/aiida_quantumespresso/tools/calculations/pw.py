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

    def get_occupations_dict(self) -> dict:
        """Return the occupations for a PwCalculation/PwBaseWorkChain node as a standard python dictionary."""

        # assert first that this is a Hubbard calculation
        try:
            type(self._node.inputs.structure) == HubbardStructureData
        except exceptions.NotExistent as exc:
            raise ValueError('The input structure is not a HubbardStructureData node.') from exc

        
        try: 
            with self._node.outputs.retrieved.base.repository.open('data-file-schema.xml') as xml_file:
                xml_parsed = ElementTree.parse(xml_file)

                schema_filepath = get_schema_filepath(xml_parsed)
                xsd = XMLSchema(schema_filepath)
                xml_dictionary = xsd.to_dict(xml_parsed, validation='skip')

                output_matrices = {}
                for  ientry, atom_dict in enumerate(xml_dictionary['output']['dft']['dftU']['Hubbard_ns']):
                    atom_specie = atom_dict['@specie']
                    shell_label = atom_dict['@label']
                    spin = 'up' if atom_dict['@spin'] == 1 else 'down'
                    shell_dims = atom_dict['@dims']
                    atom_index = atom_dict['@index']
                    occ_matrix = np.array(atom_dict['$']).reshape(shell_dims)

                    atom_label = f'Atom_{atom_index}'

                    if atom_label not in output_matrices:
                        output_matrices[atom_label] = {}
                        output_matrices[atom_label]['specie'] = atom_specie
                        output_matrices[atom_label]['shell'] = shell_label
                        output_matrices[atom_label]['occupation_matrix'] = {'up': {}, 'down': {}}
                    
                    output_matrices[atom_label]['occupation_matrix'][spin] = occ_matrix
        except Exception as exc:
            # raise just a warning and return None
            print(f'Warning: could not parse occupation matrices from XML file: {exc}')
            return None

        return output_matrices
