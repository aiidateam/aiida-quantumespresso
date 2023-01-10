# -*- coding: utf-8 -*-
"""CalcFunction to analyse the symmetry properties of a crystal structure and its atomic sites.

Returns a supercell with a marked absorbing atom for each symmetrically non-equivalent site in the system.
"""
from aiida import orm
from aiida.common import ValidationError
from aiida.common.constants import elements
from aiida.engine import calcfunction
from aiida.tools import spglib_tuple_to_structure, structure_to_spglib_tuple
import numpy as np
import spglib


@calcfunction
def get_xspectra_structures(structure, spglib_settings=None, **kwargs):  # pylint: disable=too-many-statements
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
        - standardize_structure: a Bool node defining if the input structure should
                                 standardized using spglib. The input structure will be
                                 standardized if no input is given, or if the crystal system
                                 is triclinic.
        - absorbing_elements_list: a List node defining the list of elements to consider when
                                   producing structures. All elements in the structure will be
                                   considered if no input is given.

    :param structure: the StructureData object to be analysed
    :param spglib_settings: an optional Dict object containing overrides for the symmetry
                            tolerance parameters used by spglib (symmprec, angle_tolerance).
    :returns: StructureData objects for the standardized crystal structure, the supercell, and
              all generated structure and associated symmetry data
    """

    from aiida.orm.nodes.data.structure import Kind, Site, StructureData

    elements_present = [kind.symbol for kind in structure.kinds]
    unwrapped_kwargs = {key: node for key, node in kwargs.items() if isinstance(node, orm.Data)}
    if 'abs_atom_marker' in unwrapped_kwargs.keys():
        abs_atom_marker = unwrapped_kwargs['abs_atom_marker'].value
        if abs_atom_marker in elements_present:
            raise ValidationError(
                f'The marker given for the absorbing atom ("{abs_atom_marker}") matches an existing Kind in the '
                f'input structure ({elements_present}.'
            )
        unwrapped_kwargs.pop('abs_atom_marker')
    else:
        abs_atom_marker = 'X'
    if 'supercell_min_parameter' in unwrapped_kwargs.keys():
        supercell_min_parameter = unwrapped_kwargs.pop('supercell_min_parameter').value
        if supercell_min_parameter < 0:
            raise ValueError(f'The requested minimum supercell parameter ({supercell_min_parameter}) is less than 0.')
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

    output_params = {}
    result = {}
    incoming_structure_size = len(structure.sites)
    incoming_structure_cell = structure.cell
    incoming_structure_params = structure.cell_lengths
    output_params['input_structure_num_sites'] = incoming_structure_size
    output_params['input_structure_cell_matrix'] = incoming_structure_cell
    output_params['input_structure_cell_lengths'] = incoming_structure_params

    incoming_structure_tuple = structure_to_spglib_tuple(structure)
    if spglib_settings:
        valid_keys = ['symprec', 'angle_tolerance']
        spglib_args = {Key: Value for Key, Value in spglib_settings.items() if Key in valid_keys}
    else:
        spglib_args = {}

    symmetry_dataset = spglib.get_symmetry_dataset(incoming_structure_tuple[0], **spglib_args)

    # if there is no symmetry to exploit, or no standardization is desired, then we just use
    # the input structure in the following steps. This is done to account for the case where
    # the user has submitted an improper crystal for calculation work and doesn't want it to
    # be changed.
    if symmetry_dataset['number'] in [1, 2] or not standardize_structure:
        standardized_structure_node = spglib_tuple_to_structure(incoming_structure_tuple[0])
        input_standardized = False
    else:  # otherwise, we proceed with the standardized structure.
        standardized_structure_tuple = spglib.standardize_cell(incoming_structure_tuple[0], **spglib_args)
        standardized_structure_node = spglib_tuple_to_structure(standardized_structure_tuple)
        # if we are standardizing the structure, then we need to update the symmetry
        # information for the standardized structure
        symmetry_dataset = spglib.get_symmetry_dataset(standardized_structure_tuple, **spglib_args)
        input_standardized = True

    equivalent_atoms_array = symmetry_dataset['equivalent_atoms']
    element_types = symmetry_dataset['std_types']

    equivalency_dict = {}
    for symmetry_value, element_type in zip(equivalent_atoms_array, element_types):
        if elements_defined:  # only process the elements given in the list
            if f'site_{symmetry_value}' in equivalency_dict:
                equivalency_dict[f'site_{symmetry_value}']['multiplicity'] += 1
            elif elements[element_type]['symbol'] not in abs_elements_list:
                pass
            else:
                equivalency_dict[f'site_{symmetry_value}'] = {
                    'symbol': elements[element_type]['symbol'],
                    'site_index': symmetry_value,
                    'multiplicity': 1
                }
        else:  # process everything in the system
            if f'site_{symmetry_value}' in equivalency_dict:
                equivalency_dict[f'site_{symmetry_value}']['multiplicity'] += 1
            else:
                equivalency_dict[f'site_{symmetry_value}'] = {
                    'symbol': elements[element_type]['symbol'],
                    'site_index': symmetry_value,
                    'multiplicity': 1
                }
    if elements_defined:
        output_params['absorbing_elements_list'] = abs_elements_list
    else:
        output_params['absorbing_elements_list'] = [Kind.symbol for Kind in structure.kinds]

    output_params['equivalent_sites_data'] = equivalency_dict
    output_params['spacegroup_number'] = symmetry_dataset['number']
    output_params['international_symbol'] = symmetry_dataset['international']

    result['standardized_structure'] = standardized_structure_node
    output_params['input_standardized'] = input_standardized
    if input_standardized:
        output_params['standardized_structure_num_sites'] = len(standardized_structure_node.sites)
        output_params['standardized_structure_cell_matrix'] = standardized_structure_node.cell
        output_params['standardized_structure_params'] = standardized_structure_node.cell_lengths

    standard_pbc = standardized_structure_node.cell_lengths

    multiples = []
    if not scale_unit_cell:
        multiples = [1, 1, 1]
    else:
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
    output_params['supercell_params'] = new_supercell.cell_lengths

    for value in equivalency_dict.values():
        target_site = value['site_index']
        marked_structure = StructureData()
        supercell_kinds = new_supercell.kinds

        marked_structure.set_cell(new_supercell.cell)
        for kind in supercell_kinds:
            marked_structure.append_kind(kind)

        index_counter = 0
        for site in new_supercell.sites:
            if index_counter == target_site:
                absorbing_kind = Kind(name=abs_atom_marker, symbols=site.kind_name)
                absorbing_site = Site(kind_name=absorbing_kind.name, position=site.position)
                marked_structure.append_kind(absorbing_kind)
                marked_structure.append_site(absorbing_site)
            else:
                new_site = Site(kind_name=site.kind_name, position=site.position)
                marked_structure.append_site(new_site)
            index_counter += 1

        result[f'site_{target_site}'] = marked_structure

    result['output_parameters'] = orm.Dict(dict=output_params)

    return result
