# -*- coding: utf-8 -*-
"""CalcFunction to analyse the symmetry properties of a crystal structure and its atomic sites.

Returns a supercell with a marked absorbing atom for each symmetrically non-equivalent site in the system.
"""
from aiida import orm
from aiida.common import ValidationError
from aiida.common.constants import elements
from aiida.engine import calcfunction
from aiida.orm.nodes.data.structure import Kind, Site, StructureData
from aiida.tools import spglib_tuple_to_structure, structure_to_spglib_tuple
import numpy as np
import spglib


@calcfunction
def get_xspectra_structures(structure, **kwargs):  # pylint: disable=too-many-statements
    """Read a StructureData object using spglib for its symmetry data and return structures for XSpectra calculations.

    Takes an incoming StructureData node and prepares structures suitable for calculation with
    xspectra.x. The basic criteria for a suitable structure is for the cell parameters to be
    sufficiently large to reduce spurious interactions between localised core-holes in
    neighbouring periodic images and for one atomic site to be marked as the absorbing atom
    site.

    The default setting for the cell size is 8.0 angstrom, which can be changed by including
    'supercell_min_parameter' as a Float node in the keyword arguments. Accepted keyword
    arguments are:

        - abs_atom_marker: a Str node defining the name of the absorbing atom Kind. The
                           absorbing Kind will be labelled 'X' if no input is given.
        - supercell_min_parameter: a Float node defining the minimum cell length in
                                   angstrom for the resulting supercell, and thus all output
                                   structures. The default value of 8.0 angstrom will be used
                                   if no input is given. Setting this value to 0.0 will
                                   instruct the CF to not scale up the input structure.
        - standardize_structure: a Bool object defining if the input structure should
                                 standardized using spglib. The input structure will be
                                 standardized if no input is given, or if the crystal system
                                 is triclinic.
        - absorbing_elements_list: a List object defining the list of elements to consider
                                   when producing structures. All elements in the structure
                                   will be considered if no input is given.
        - is_molecule_input: a Bool object to define for the function that the input structure is
                    a molecule and not a periodic solid system. Required in order to instruct
                    the CF to use Pymatgen rather than spglib to determine the symmetry. The
                    CF will assume the structure to be a periodic solid if no input is given.
        - spglib_settings: an optional Dict object containing overrides for the symmetry
                            tolerance parameters used by spglib (symmprec, angle_tolerance).
        - pymatgen_settings: an optional Dict object containing overrides for the symmetry
                             tolerance parameters used by Pymatgen when processing molecular
                             systems (tolerance, eigen_tolerance, matrix_tolerance).

    :param structure: the StructureData object to be analysed
    :returns: StructureData objects for the standardized crystal structure, the supercell, and
              all generated structure and associated symmetry data
    """

    is_molecule_input = False
    elements_present = [kind.symbol for kind in structure.kinds]
    unwrapped_kwargs = {key: node for key, node in kwargs.items() if isinstance(node, orm.Data)}
    if 'abs_atom_marker' in unwrapped_kwargs.keys():
        abs_atom_marker = unwrapped_kwargs['abs_atom_marker'].value
        if abs_atom_marker in elements_present:
            raise ValidationError(
                f'The marker given for the absorbing atom ("{abs_atom_marker}") should not match an existing Kind in '
                f'the input structure ({elements_present}.'
            )
        unwrapped_kwargs.pop('abs_atom_marker')
    else:
        abs_atom_marker = 'X'
    if 'supercell_min_parameter' in unwrapped_kwargs.keys():
        supercell_min_parameter = unwrapped_kwargs.pop('supercell_min_parameter').value
        if supercell_min_parameter < 0:
            raise ValueError(
                f'The requested minimum supercell parameter ({supercell_min_parameter}) should not be'
                ' less than 0.'
            )
        elif supercell_min_parameter == 0:  # useful if no core-hole treatment is required
            scale_unit_cell = False
        else:
            scale_unit_cell = True
    else:
        supercell_min_parameter = 8.0  # angstroms
        scale_unit_cell = True
    if 'standardize_structure' in unwrapped_kwargs.keys():
        standardize_structure = unwrapped_kwargs['standardize_structure'].value
        unwrapped_kwargs.pop('standardize_structure')
    else:
        standardize_structure = True
    if 'absorbing_elements_list' in unwrapped_kwargs.keys():
        elements_defined = True
        abs_elements_list = unwrapped_kwargs['absorbing_elements_list'].get_list()
        # confirm that the elements requested are actually in the input structure
        for req_element in abs_elements_list:
            if req_element not in elements_present:
                raise ValidationError(
                    f'Requested elements {abs_elements_list} do not match those in the given structure:'
                    f' {elements_present}'
                )
    else:
        elements_defined = False
        abs_elements_list = [Kind.symbol for Kind in structure.kinds]
    if 'is_molecule_input' in unwrapped_kwargs.keys():
        is_molecule_input = unwrapped_kwargs['is_molecule_input'].value
        # If we are working with a molecule, check for pymatgen_settings
        if 'pymatgen_settings' in unwrapped_kwargs.keys():
            pymatgen_settings_dict = unwrapped_kwargs.keys()['pymatgen_settings'].get_dict()
            valid_keys = ['tolerance', 'eigen_tolerance', 'matrix_tolerance']
            pymatgen_kwargs = {key: value for key, value in pymatgen_settings_dict.items() if key in valid_keys}
        else:
            pymatgen_kwargs = {}

    if 'spglib_settings' in unwrapped_kwargs.keys():
        spglib_settings_dict = unwrapped_kwargs['spglib_settings'].get_dict()
        valid_keys = ['symprec', 'angle_tolerance']
        spglib_kwargs = {key: value for key, value in spglib_settings_dict.items() if key in valid_keys}
    else:
        spglib_kwargs = {}

    output_params = {}
    result = {}

    output_params['absorbing_elements_list'] = abs_elements_list

    incoming_structure_size = len(structure.sites)
    incoming_structure_cell = structure.cell
    incoming_structure_params = structure.cell_lengths

    output_params['input_structure_num_sites'] = incoming_structure_size
    output_params['input_structure_cell_matrix'] = incoming_structure_cell
    output_params['input_structure_cell_lengths'] = incoming_structure_params

    # Process a non-periodic system
    if is_molecule_input:
        from pymatgen.symmetry.analyzer import PointGroupAnalyzer

        pymatgen_molecule = structure.get_pymatgen_molecule()
        centered_molecule = pymatgen_molecule.get_centered_molecule()
        centered_positions = centered_molecule.cart_coords

        analyzer_data = PointGroupAnalyzer(pymatgen_molecule, **pymatgen_kwargs)
        eq_atoms_data = analyzer_data.get_equivalent_atoms()['eq_sets']
        atomic_nos = pymatgen_molecule.atomic_numbers
        point_group = str(analyzer_data.get_pointgroup())
        output_params['molecule_point_group'] = point_group

        equivalency_dict = {}

        for key in eq_atoms_data:
            site_symbol = elements[atomic_nos[key]]['symbol']
            if site_symbol in abs_elements_list:
                equivalency_dict[f'site_{key}'] = {}
                atom_no_set = eq_atoms_data[key]
                equivalency_dict[f'site_{key}']['site_index'] = key
                equivalency_dict[f'site_{key}']['equivalent_sites_list'] = atom_no_set
                equivalency_dict[f'site_{key}']['multiplicity'] = len(atom_no_set)
                equivalency_dict[f'site_{key}']['symbol'] = site_symbol

        output_params['equivalent_sites_data'] = equivalency_dict

        x_pos = []
        y_pos = []
        z_pos = []
        for site in structure.sites:
            x_pos.append(site.position[0])
            y_pos.append(site.position[1])
            z_pos.append(site.position[2])

        # We assume that the Martyna-Tuckerman approach will be used in the core-hole SCF
        # thus, at a minimum, the new box for the molecule needs to be 2x the extent of the
        # molecule in x, y, and z. We therefore use the larger of either the MT distance or
        # `supercell_min_parameter` for each cell parameter.
        high_x = x_pos[(np.argmax(x_pos))]
        high_y = y_pos[(np.argmax(y_pos))]
        high_z = z_pos[(np.argmax(z_pos))]

        low_x = x_pos[(np.argmin(x_pos))]
        low_y = y_pos[(np.argmin(y_pos))]
        low_z = z_pos[(np.argmin(z_pos))]

        box_a = max((high_x - low_x) * 2, supercell_min_parameter)
        box_b = max((high_y - low_y) * 2, supercell_min_parameter)
        box_c = max((high_z - low_z) * 2, supercell_min_parameter)

        new_supercell = StructureData()
        new_supercell.set_cell(value=([box_a, 0., 0.], [0., box_b, 0.], [0., 0., box_c]))

        for kind in structure.kinds:
            new_supercell.append_kind(kind)

        for site in structure.sites:
            new_supercell.append_site(site)

        # Reset the positions so that the centre of mass is at the centre of the new box
        new_supercell.reset_sites_positions(centered_positions + (np.array(new_supercell.cell_lengths) / 2))

        output_params['supercell_cell_matrix'] = new_supercell.cell
        output_params['supercell_cell_lengths'] = new_supercell.cell_lengths
        output_params['supercell_num_sites'] = len(new_supercell.sites)
        result['supercell'] = new_supercell

    # Process a periodic system
    else:
        incoming_structure_tuple = structure_to_spglib_tuple(structure)

        symmetry_dataset = spglib.get_symmetry_dataset(incoming_structure_tuple[0], **spglib_kwargs)

        # if there is no symmetry to exploit, or no standardization is desired, then we just use
        # the input structure in the following steps. This is done to account for the case where
        # the user has submitted an improper crystal for calculation work and doesn't want it to
        # be changed.
        if symmetry_dataset['number'] in [1, 2] or not standardize_structure:
            standardized_structure_node = spglib_tuple_to_structure(incoming_structure_tuple[0])
            structure_is_standardized = False
        else:  # otherwise, we proceed with the standardized structure.
            standardized_structure_tuple = spglib.standardize_cell(incoming_structure_tuple[0], **spglib_kwargs)
            standardized_structure_node = spglib_tuple_to_structure(standardized_structure_tuple)
            # if we are standardizing the structure, then we need to update the symmetry
            # information for the standardized structure
            symmetry_dataset = spglib.get_symmetry_dataset(standardized_structure_tuple, **spglib_kwargs)
            structure_is_standardized = True

        equivalent_atoms_array = symmetry_dataset['equivalent_atoms']
        element_types = symmetry_dataset['std_types']

        equivalency_dict = {}
        index_counter = 0
        for symmetry_value, element_type in zip(equivalent_atoms_array, element_types):
            if elements_defined:  # only process the elements given in the list
                if f'site_{symmetry_value}' in equivalency_dict:
                    equivalency_dict[f'site_{symmetry_value}']['equivalent_sites_list'].append(index_counter)
                elif elements[element_type]['symbol'] not in abs_elements_list:
                    pass
                else:
                    equivalency_dict[f'site_{symmetry_value}'] = {
                        'symbol': elements[element_type]['symbol'],
                        'site_index': symmetry_value,
                        'equivalent_sites_list': [symmetry_value]
                    }
            else:  # process everything in the system
                if f'site_{symmetry_value}' in equivalency_dict:
                    equivalency_dict[f'site_{symmetry_value}']['equivalent_sites_list'].append(index_counter)
                else:
                    equivalency_dict[f'site_{symmetry_value}'] = {
                        'symbol': elements[element_type]['symbol'],
                        'site_index': symmetry_value,
                        'equivalent_sites_list': [symmetry_value]
                    }
            index_counter += 1

        for value in equivalency_dict.values():
            value['multiplicity'] = len(value['equivalent_sites_list'])

        output_params['equivalent_sites_data'] = equivalency_dict
        output_params['spacegroup_number'] = symmetry_dataset['number']
        output_params['international_symbol'] = symmetry_dataset['international']

        output_params['structure_is_standardized'] = structure_is_standardized
        if structure_is_standardized:
            result['standardized_structure'] = standardized_structure_node
            output_params['standardized_structure_num_sites'] = len(standardized_structure_node.sites)
            output_params['standardized_structure_cell_matrix'] = standardized_structure_node.cell
            output_params['standardized_structure_params'] = standardized_structure_node.cell_lengths

        standard_pbc = standardized_structure_node.cell_lengths

        multiples = [1, 1, 1]
        if scale_unit_cell:
            multiples = []
            for length in standard_pbc:
                multiple = np.ceil(supercell_min_parameter / length)
                multiples.append(int(multiple))

        ase_structure = standardized_structure_node.get_ase()
        ase_supercell = ase_structure * multiples
        new_supercell = StructureData(ase=ase_supercell)

        result['supercell'] = new_supercell
        output_params['supercell_factors'] = multiples
        output_params['supercell_num_sites'] = len(new_supercell.sites)
        output_params['supercell_cell_matrix'] = new_supercell.cell
        output_params['supercell_cell_lengths'] = new_supercell.cell_lengths

    for value in equivalency_dict.values():
        target_site = value['site_index']
        marked_structure = StructureData()
        supercell_kinds = {kind.name: kind for kind in new_supercell.kinds}
        marked_structure.set_cell(new_supercell.cell)

        for index, site in enumerate(new_supercell.sites):
            if index == target_site:
                absorbing_kind = Kind(name=abs_atom_marker, symbols=site.kind_name)
                absorbing_site = Site(kind_name=absorbing_kind.name, position=site.position)
                marked_structure.append_kind(absorbing_kind)
                marked_structure.append_site(absorbing_site)
            else:
                if site.kind_name not in [kind.name for kind in marked_structure.kinds]:
                    marked_structure.append_kind(supercell_kinds[site.kind_name])
                new_site = Site(kind_name=site.kind_name, position=site.position)
                marked_structure.append_site(new_site)

        result[f'site_{target_site}_{value["symbol"]}'] = marked_structure

    output_params['is_molecule_input'] = is_molecule_input
    result['output_parameters'] = orm.Dict(dict=output_params)

    return result
