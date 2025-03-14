# -*- coding: utf-8 -*-
"""CalcFunction to analyse the symmetry properties of a crystal structure and its atomic sites.

Returns a supercell with a marked absorbing atom for each symmetrically non-equivalent site in the system.
"""
import warnings

from aiida import orm
from aiida.common import ValidationError
from aiida.common.constants import elements
from aiida.engine import calcfunction
from aiida.orm.nodes.data.structure import Kind, Site, StructureData
from aiida.tools import spglib_tuple_to_structure, structure_to_spglib_tuple
import numpy as np
import spglib

from aiida_quantumespresso.utils.hubbard import HubbardStructureData, HubbardUtils

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


def process_molecule_input(structure, **kwargs):  # pylint: disable=too-many-statements
    """Generate symmetry data and supercell suitable for molecular systems using Pymatgen."""
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer

    result = {'output_params': {}}
    if 'pymatgen_settings' in kwargs:
        pymatgen_settings_dict = kwargs['pymatgen_settings'].get_dict()
        valid_keys = ['tolerance', 'eigen_tolerance', 'matrix_tolerance']
        pymatgen_kwargs = {key: value for key, value in pymatgen_settings_dict.items() if key in valid_keys}
    else:
        pymatgen_kwargs = {}

    abs_elements_list = kwargs.get('abs_elements_list', [kind.symbol for kind in structure.kinds])
    pymatgen_molecule = structure.get_pymatgen_molecule()
    centered_molecule = pymatgen_molecule.get_centered_molecule()
    centered_positions = centered_molecule.cart_coords

    analyzer_data = PointGroupAnalyzer(pymatgen_molecule, **pymatgen_kwargs)
    eq_atoms_data = analyzer_data.get_equivalent_atoms()['eq_sets']
    atomic_nos = pymatgen_molecule.atomic_numbers
    point_group = str(analyzer_data.get_pointgroup())
    result['output_params']['molecule_point_group'] = point_group

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

    result['output_params']['equivalent_sites_data'] = equivalency_dict

    supercell_min_parameter = kwargs.get('supercell_min_parameter', 8.0)

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

    result['output_params']['supercell_cell_matrix'] = new_supercell.cell
    result['output_params']['supercell_cell_lengths'] = new_supercell.cell_lengths
    result['output_params']['supercell_num_sites'] = len(new_supercell.sites)
    result['supercell'] = new_supercell

    return result


def get_spglib_equivalency_dict(equivalent_atoms_array, abs_elements_list, element_types, type_mapping_dict):
    """Convert the equivalent atoms array into the required dictionary format."""

    equivalency_dict = {}
    for index, symmetry_types in enumerate(zip(equivalent_atoms_array, element_types)):
        symmetry_value, element_type = symmetry_types
        if type_mapping_dict[str(element_type)].symbol in abs_elements_list:
            if f'site_{symmetry_value}' in equivalency_dict:
                equivalency_dict[f'site_{symmetry_value}']['equivalent_sites_list'].append(index)
            else:
                equivalency_dict[f'site_{symmetry_value}'] = {
                    'kind_name': type_mapping_dict[str(element_type)].name,
                    'symbol': type_mapping_dict[str(element_type)].symbol,
                    'site_index': symmetry_value,
                    'equivalent_sites_list': [symmetry_value]
                }

    return equivalency_dict


def get_supercell(structure, supercell_min_parameter, is_hubbard_structure, **kwargs):
    """Generate a supercell of a PBC system, given a minimum cell length requirement."""
    standard_pbc = structure.cell_lengths

    result = {}
    multiples = []
    for length in standard_pbc:
        multiple = np.ceil(supercell_min_parameter / length)
        multiples.append(int(multiple))

    result['multiples'] = multiples
    ase_structure = structure.get_ase()
    ase_supercell = ase_structure * multiples

    # Here we construct the supercell site-by-site in order to keep the
    # correct Kind ordering (if any)
    blank_supercell = StructureData(ase=ase_supercell)
    new_supercell = StructureData()
    new_supercell.set_cell(blank_supercell.cell)
    num_extensions = np.product(multiples)
    types_order = kwargs['types_order']
    type_mapping_dict = kwargs['type_mapping_dict']
    supercell_types_order = []

    for _ in range(0, num_extensions):
        for type_number in types_order:
            supercell_types_order.append(type_number)

    for site, type_number in zip(blank_supercell.sites, supercell_types_order):
        kind_present = type_mapping_dict[str(type_number)]
        if kind_present.name not in [kind.name for kind in new_supercell.kinds]:
            new_supercell.append_kind(kind_present)
        new_site = Site(kind_name=kind_present.name, position=site.position)
        new_supercell.append_site(new_site)

    if is_hubbard_structure:
        # Scale up the hubbard parameters to match and return the HubbardStructureData.
        # We can exploit the fact that `get_hubbard_for_supercell` will return a HubbardStructureData node
        # with the same hubbard parameters in the case where the input structure and the supercell are the
        # same (i.e. multiples == [1, 1, 1])
        hubbard_manip = HubbardUtils(structure)
        new_hubbard_supercell = hubbard_manip.get_hubbard_for_supercell(new_supercell)
        new_supercell = new_hubbard_supercell
        supercell_hubbard_params = new_supercell.hubbard
        result['supercell_hubbard_params'] = supercell_hubbard_params

    result['new_supercell'] = new_supercell
    return result


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
        - use_element_types: a Bool object to indicate that symmetry analysis should consider all
                             sites of the same element to be equal and ignore any special Kind names
                             from the parent structure. For instance, use_element_types = True would
                             consider sites for Kinds 'Si' and 'Si1' to be equivalent if both are sites
                             containing silicon. Defaults to True.
        - parse_symmetry: a Bool object to indicate that the symmetry data should be generated using spglib
                          for the input structure. If False, will instead use symmetry data provided
                          via the `equivalent_sites_data` and `spacegroup_number` inputs and generate all
                          the structures based on this information. Defaults to True.
        - equivalent_sites_data: a Dict object taking the form of the `equivalent_sites_data` dict normally
                                 normally generated by the CF. Only read if `parse_symmetry = False`.
                                 Example: {'site_0' : {'symbol' : 'Si', 'site_index' : 0, 'multiplicity' : 8}}
                                 Note that, for each site, only 'symbol', 'site_index', and 'multiplicity' are
                                 required, all other options are ignored.
        - spacegroup_number: an Int object defining the spacegroup number of the input structure, only read
                             if `parse_symmetry = False`.
        - spglib_settings: an optional Dict object containing overrides for the symmetry
                            tolerance parameters used by spglib (symmprec, angle_tolerance).
        - pymatgen_settings: an optional Dict object containing overrides for the symmetry
                             tolerance parameters used by Pymatgen when processing molecular
                             systems (tolerance, eigen_tolerance, matrix_tolerance).

    :param structure: the StructureData object to be analysed
    :returns: StructureData objects for the standardized crystal structure, the supercell, and
              all generated structure and associated symmetry data
    """

    elements_present = [kind.symbol for kind in structure.kinds]
    unwrapped_kwargs = {key: node for key, node in kwargs.items() if isinstance(node, orm.Data)}

    is_molecule_input = unwrapped_kwargs.get('is_molecule_input', False)
    use_element_types = unwrapped_kwargs.get('use_element_types', True)
    parse_symmetry = unwrapped_kwargs.get('parse_symmetry', True)

    if 'abs_atom_marker' in unwrapped_kwargs:
        abs_atom_marker = unwrapped_kwargs['abs_atom_marker'].value
        if abs_atom_marker in elements_present:
            raise ValidationError(
                f'The marker given for the absorbing atom ("{abs_atom_marker}") should not match an existing Kind in '
                f'the input structure ({elements_present}.'
            )
        unwrapped_kwargs.pop('abs_atom_marker')
    else:
        abs_atom_marker = 'X'
    if 'supercell_min_parameter' in unwrapped_kwargs:
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
    if 'standardize_structure' in unwrapped_kwargs:
        standardize_structure = unwrapped_kwargs['standardize_structure'].value
        unwrapped_kwargs.pop('standardize_structure')
    else:
        standardize_structure = True
    if 'absorbing_elements_list' in unwrapped_kwargs:
        abs_elements_list = unwrapped_kwargs['absorbing_elements_list'].get_list()
        # confirm that the elements requested are actually in the input structure
        for req_element in abs_elements_list:
            if req_element not in elements_present:
                raise ValidationError(
                    f'Requested elements {abs_elements_list} do not match those in the given structure:'
                    f' {elements_present}'
                )
    else:
        abs_elements_list = []
        for kind in structure.kinds:
            if kind.symbol not in abs_elements_list:
                abs_elements_list.append(kind.symbol)

    if 'spglib_settings' in unwrapped_kwargs:
        spglib_settings_dict = unwrapped_kwargs['spglib_settings'].get_dict()
        valid_keys = ['symprec', 'angle_tolerance']
        spglib_kwargs = {key: value for key, value in spglib_settings_dict.items() if key in valid_keys}
    else:
        spglib_kwargs = {}

    if isinstance(structure, HubbardStructureData):
        is_hubbard_structure = True
        if standardize_structure:
            raise ValidationError(
                'Incoming structure set to be standardized, but hubbard data has been found. '
                'Please set ``standardize_structure`` to false in ``**kwargs`` to preserve the hubbard data.'
            )
    else:
        is_hubbard_structure = False

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
    get_supercell_inputs = {'is_hubbard_structure': is_hubbard_structure}

    incoming_structure_tuple = structure_to_spglib_tuple(structure)
    spglib_tuple = incoming_structure_tuple[0]
    types_order = spglib_tuple[-1]
    kinds_information = incoming_structure_tuple[1]
    kinds_list = incoming_structure_tuple[2]

    # We need a way to reliably convert type number into element, so we
    # first create a mapping of assigned number to kind name then a mapping
    # of assigned number to ``Kind``
    # Note that we set `types_order` and `types_mapping_dict` using the
    # input structure in case the symmetry is not parsed later. If it is,
    # then `types_order` and `type_mapping_dict` will be updated as needed.

    type_name_mapping = {str(value): key for key, value in kinds_information.items()}
    type_mapping_dict = {}

    for key, value in type_name_mapping.items():
        for kind in kinds_list:
            if value == kind.name:
                type_mapping_dict[key] = kind

    get_supercell_inputs['types_order'] = types_order
    get_supercell_inputs['type_mapping_dict'] = type_mapping_dict

    if is_molecule_input:
        process_molecule_result = process_molecule_input(structure, **kwargs)
        new_supercell = process_molecule_result['supercell']
        for key, value in process_molecule_result['output_params'].items():
            output_params[key] = value
        equivalency_dict = output_params['equivalent_sites_data']
    # Process a periodic system
    elif not parse_symmetry:
        equivalency_dict = kwargs['equivalent_sites_data'].get_dict()
        output_params['equivalent_sites_data'] = equivalency_dict
        output_params['spacegroup_number'] = kwargs['spacegroup_number'].value
        output_params['symmetry_parsed'] = False
        output_params['structure_is_standardized'] = False
        get_supercell_inputs['structure'] = structure
        get_supercell_inputs['supercell_min_parameter'] = supercell_min_parameter

        if not scale_unit_cell:
            multiples = [1, 1, 1]
            new_supercell = structure
        else:
            get_supercell_result = get_supercell(**get_supercell_inputs)
            multiples = get_supercell_result['multiples']
            new_supercell = get_supercell_result['new_supercell']
            output_params['supercell_factors'] = multiples

        result['supercell'] = new_supercell
        output_params['supercell_num_sites'] = len(new_supercell.sites)
        output_params['supercell_cell_matrix'] = new_supercell.cell
        output_params['supercell_cell_lengths'] = new_supercell.cell_lengths

    else:
        # By default, `structure_to_spglib_tuple` gives different
        # ``Kinds`` of the same element a distinct atomic number by
        # multiplying the normal atomic number by 1000, then adding
        # 100 for each distinct duplicate. if we want to treat all sites
        # of the same element as equal, then we must therefore briefly
        # operate on a "cleaned" version of the structure tuple where this
        # new label is reduced to its normal element number.
        if use_element_types:
            cleaned_structure_tuple = (spglib_tuple[0], spglib_tuple[1], [])
            for i in spglib_tuple[2]:
                if i >= 1000:
                    new_i = int(np.trunc(i / 1000))
                else:
                    new_i = i
                cleaned_structure_tuple[2].append(new_i)
            symmetry_dataset = spglib.get_symmetry_dataset(cleaned_structure_tuple, **spglib_kwargs)
        else:
            symmetry_dataset = spglib.get_symmetry_dataset(spglib_tuple, **spglib_kwargs)

        # if there is no symmetry to exploit, or no standardization is desired, then we just use
        # the input structure in the following steps. This is done to account for the case where
        # the user has submitted an improper crystal for calculation work and doesn't want it to
        # be changed.
        if symmetry_dataset['number'] in [1, 2] or not standardize_structure:
            standardized_structure_node = structure
            structure_is_standardized = False
        else:  # otherwise, we proceed with the standardized structure.
            standardized_structure_tuple = spglib.standardize_cell(spglib_tuple, **spglib_kwargs)
            standardized_structure_node = spglib_tuple_to_structure(
                standardized_structure_tuple, kinds_information, kinds_list
            )
            # if we are standardizing the structure, then we need to update the symmetry
            # information for the standardized structure
            symmetry_dataset = spglib.get_symmetry_dataset(standardized_structure_tuple, **spglib_kwargs)
            structure_is_standardized = True

        get_supercell_inputs['structure'] = standardized_structure_node
        get_supercell_inputs['supercell_min_parameter'] = supercell_min_parameter
        equivalent_atoms_array = symmetry_dataset['equivalent_atoms']

        if structure_is_standardized:
            element_types = symmetry_dataset['std_types']
        else:  # convert the 'std_types' from standardized to primitive cell
            # we generate the type-specific data on-the-fly since we need to
            # know which type (and thus kind) *should* be at each site
            # even if we "cleaned" the structure previously
            non_cleaned_dataset = spglib.get_symmetry_dataset(spglib_tuple, **spglib_kwargs)
            spglib_std_types = non_cleaned_dataset['std_types']
            spglib_map_to_prim = non_cleaned_dataset['mapping_to_primitive']
            spglib_std_map_to_prim = non_cleaned_dataset['std_mapping_to_primitive']

            map_std_pos_to_type = {}
            for position, atom_type in zip(spglib_std_map_to_prim, spglib_std_types):
                map_std_pos_to_type[str(position)] = atom_type
            primitive_types = []
            for position in spglib_map_to_prim:
                atom_type = map_std_pos_to_type[str(position)]
                primitive_types.append(atom_type)
            element_types = primitive_types

        equivalency_dict = get_spglib_equivalency_dict(
            equivalent_atoms_array, abs_elements_list, element_types, type_mapping_dict
        )
        output_params['symmetry_parsed'] = True
        # update the types order, in case it has changed during symmetry analysis
        get_supercell_inputs['types_order'] = element_types

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

        if not scale_unit_cell:
            multiples = [1, 1, 1]
            new_supercell = standardized_structure_node
        else:
            get_supercell_result = get_supercell(**get_supercell_inputs)
            multiples = get_supercell_result['multiples']
            new_supercell = get_supercell_result['new_supercell']

        result['supercell'] = new_supercell

        output_params['supercell_factors'] = multiples
        output_params['supercell_num_sites'] = len(new_supercell.sites)
        output_params['supercell_cell_matrix'] = new_supercell.cell
        output_params['supercell_cell_lengths'] = new_supercell.cell_lengths

    for value in equivalency_dict.values():
        target_site = value['site_index']
        marked_structure = StructureData()
        marked_structure.set_cell(new_supercell.cell)
        new_kind_names = [kind.name for kind in new_supercell.kinds]

        for index, site in enumerate(new_supercell.sites):
            kind_at_position = new_supercell.kinds[new_kind_names.index(site.kind_name)]
            if index == target_site:
                absorbing_kind = Kind(name=abs_atom_marker, symbols=kind_at_position.symbol)
                absorbing_site = Site(kind_name=absorbing_kind.name, position=site.position)
                marked_structure.append_kind(absorbing_kind)
                marked_structure.append_site(absorbing_site)
            else:
                if kind_at_position.name not in [kind.name for kind in marked_structure.kinds]:
                    marked_structure.append_kind(kind_at_position)
                new_site = Site(kind_name=site.kind_name, position=site.position)
                marked_structure.append_site(new_site)

        if is_hubbard_structure:
            marked_hubbard_structure = HubbardStructureData.from_structure(marked_structure)
            marked_hubbard_structure.hubbard = get_supercell_result['supercell_hubbard_params']
            result[f'site_{target_site}_{value["symbol"]}'] = marked_hubbard_structure
        else:
            result[f'site_{target_site}_{value["symbol"]}'] = marked_structure

    output_params['is_molecule_input'] = is_molecule_input
    result['output_parameters'] = orm.Dict(dict=output_params)

    return result
