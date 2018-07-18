# -*- coding: utf-8 -*-
import xmlschema
from defusedxml import ElementTree

from aiida_quantumespresso.parsers.constants import ry_to_ev, hartree_to_ev, bohr_to_ang, ry_si, bohr_si
from .versions import get_schema_filepath, get_default_schema_filepath


def cell_volume(a1, a2, a3):
    r"""Returns the volume of the primitive cell: :math:`|\vec a_1\cdot(\vec a_2\cross \vec a_3)|`"""
    a_mid_0 = a2[1] * a3[2] - a2[2] * a3[1]
    a_mid_1 = a2[2] * a3[0] - a2[0] * a3[2]
    a_mid_2 = a2[0] * a3[1] - a2[1] * a3[0]

    return abs(float(a1[0] * a_mid_0 + a1[1] * a_mid_1 + a1[2] * a_mid_2))


def parse_pw_xml_post_6_2(xml_file):
    """
    """
    try:
        xml = ElementTree.parse(xml_file)
    except IOError:
        raise ValueError('could not open and or parse the XML file {}'.format(xml_file))

    schema_filepath = get_schema_filepath(xml)

    try:
        xsd = xmlschema.XMLSchema(schema_filepath)
    except URLError:

        # If loading the XSD file specified in the XML file fails, we try the default
        schema_filepath = get_default_schema_filepath()

        try:
            xsd = xmlschema.XMLSchema(schema_filepath)
        except URLError:
            raise ValueError('could not open and or parse the XSD file {}'.format(schema_filepath))

    result = xsd.to_dict(xml, validation='lax')

    if isinstance(result, tuple):
        xml_dictionary, errors = result
    else:
        xml_dictionary = result
        errors = []

    import json
    print json.dumps(xml_dictionary, indent=4)

    lattice_vectors = [
        map(lambda x: x * bohr_to_ang, xml_dictionary['output']['atomic_structure']['cell']['a1']),
        map(lambda x: x * bohr_to_ang, xml_dictionary['output']['atomic_structure']['cell']['a2']),
        map(lambda x: x * bohr_to_ang, xml_dictionary['output']['atomic_structure']['cell']['a3']),
    ]

    if ('electric_field' in xml_dictionary['input'] and
        'electric_potential' in xml_dictionary['input']['electric_field'] and
        xml_dictionary['input']['electric_field']['electric_potential'] == 'sawtooth_potential'):
        has_electric_field = True
    else:
        has_electric_field = False

    if ('electric_field' in xml_dictionary['input'] and
        'dipole_correction' in xml_dictionary['input']['electric_field']):
        has_dipole_correction = xml_dictionary['input']['electric_field']['dipole_correction']
    else:
        has_dipole_correction = False

    if 'bands' in xml_dictionary['input'] and 'occupations' in xml_dictionary['input']['bands']:

        occupations = xml_dictionary['input']['bands']['occupations']

        if occupations == 'from_input':
            fixed_occupations = True
        else:
            fixed_occupations = False

        if 'tetrahedra' in occupations:
            tetrahedron_method = True
        else:
            tetrahedron_method = False

        if occupations == 'from_input':
            smearing_method = True
        else:
            smearing_method = False

    starting_magnetization = []
    magnetization_angle1 = []
    magnetization_angle2 = []

    for specie in xml_dictionary['output']['atomic_species']['species']:
        starting_magnetization.append(specie.get('starting_magnetization', 0.0))
        magnetization_angle1.append(specie.get('magnetization_angle1', 0.0))
        magnetization_angle2.append(specie.get('magnetization_angle2', 0.0))

    constraint_mag = 0
    if ('spin_constraints' in xml_dictionary['input'] and
        'spin_constraints' in xml_dictionary['input']['spin_constraints']):
        spin_constraints = xml_dictionary['input']['spin_constraints']['spin_constraints']

        if spin_constraints == 'atomic':
            constraint_mag = 1
        elif spin_constraints == 'atomic direction':
            constraint_mag = 2
        elif spin_constraints == 'total':
            constraint_mag = 3
        elif spin_constraints == 'total direction':
            constraint_mag = 6

    lsda = xml_dictionary['input']['spin']['lsda']
    spin_orbit_calculation = xml_dictionary['input']['spin']['spinorbit']
    non_colinear_calculation = xml_dictionary['output']['magnetization']['noncolin']
    do_magnetization = xml_dictionary['output']['magnetization']['do_magnetization']

    if non_colinear_calculation and do_magnetization:
        time_reversal = False
    else:
        time_reversal = True

    # If no specific tags are present, the default is 1
    if non_colinear_calculation or spin_orbit_calculation:
        nspin = 4
    elif lsda:
        nspin = 2
    else:
        nspin = 1

    # Detect presence of inversion symmetry, which is the case if a `crystal` symmetry has the attribute `inversion`
    # Note that `lattice` symmetries will always have inversion symmetry and should therefore be ignored
    symmetries = []
    inversion_symmetry = False

    for symmetry in xml_dictionary['output']['symmetries']['symmetry']:

        if symmetry['info']['$'] == 'crystal_symmetry' and 'inv' in symmetry['info']['@name']:
            inversion_symmetry = True
            break

        sym = {
            'rotation': [
                symmetry['rotation']['$'][0:3],
                symmetry['rotation']['$'][3:6],
                symmetry['rotation']['$'][6:9],
            ],
            'name': symmetry['info']['@name'],
            't_rev': '0'
        }

        try:
            sym['equivalent_atoms'] = symmetry['equivalent_atoms']['$']
        except KeyError:
            pass

        try:
            sym['fractional_translation'] = symmetry['fractional_translation']
        except KeyError:
            pass

        symmetries.append(sym)

    xml_data = {
        # 'pp_check_flag': true,  Currently not printed in the schema
        # 'beta_real_space': false,  Currently not printed in the schema
        'xml_warnings': [],
        'rho_cutoff_units': 'eV',
        'wfc_cutoff_units': 'eV',
        'fermi_energy_units': 'eV',
        'k_points_units': '2 pi / Angstrom',
        'symmetries_units': 'crystal',
        'constraint_mag': constraint_mag,
        'fixed_occupations': fixed_occupations,
        'tetrahedron_method': tetrahedron_method,
        'smearing_method': smearing_method,
        'magnetization_angle2': magnetization_angle2,
        'magnetization_angle1': magnetization_angle1,
        'starting_magnetization': starting_magnetization,
        'has_electric_field': has_electric_field,
        'has_dipole_correction': has_dipole_correction,
        'lda_plus_u_calculation': 'dftU' in xml_dictionary['output'],
        'format_version': xml_dictionary['general_info']['xml_format']['@NAME'],
        'format_name': xml_dictionary['general_info']['xml_format']['@VERSION'],
        'monkhorst_pack_offset': xml_dictionary['input']['k_points_IBZ']['monkhorst_pack'].values()[1:4],
        'monkhorst_pack_grid': xml_dictionary['input']['k_points_IBZ']['monkhorst_pack'].values()[4:7],
        'non_colinear_calculation': non_colinear_calculation,
        'do_magnetization': do_magnetization,
        'time_reversal_flag': time_reversal,
        'symmetries': symmetries,
        'k_points': [eigenvalue['k_point']['$'] for eigenvalue in xml_dictionary['output']['band_structure']['ks_energies']],
        'do_not_use_time_reversal': xml_dictionary['input']['symmetry_flags']['noinv'],
        'rho_cutoff': xml_dictionary['input']['basis']['ecutrho'] * hartree_to_ev,
        'number_of_bands': xml_dictionary['output']['band_structure']['nbnd'],
        'spin_orbit_domag': xml_dictionary['output']['magnetization']['do_magnetization'],
        'fft_grid': xml_dictionary['output']['basis_set']['fft_grid'].values(),
        'lsda': lsda,
        'number_of_spin_components': nspin,
        'no_time_rev_operations': xml_dictionary['input']['symmetry_flags']['no_t_rev'],
        'inversion_symmetry': inversion_symmetry,
        'number_of_bravais_symmetries': xml_dictionary['output']['symmetries']['nrot'],
        'number_of_symmetries': xml_dictionary['output']['symmetries']['nsym'],
        'wfc_cutoff': xml_dictionary['input']['basis']['ecutwfc'] * hartree_to_ev,
        'smooth_fft_grid': xml_dictionary['output']['basis_set']['fft_smooth'].values(),
        'dft_exchange_correlation': xml_dictionary['input']['dft']['functional'],
        'spin_orbit_calculation': spin_orbit_calculation,
        'creator_name': xml_dictionary['general_info']['creator']['@NAME'].lower(),
        'number_of_atomic_wfc': xml_dictionary['output']['band_structure']['num_of_atomic_wfc'],
        'creator_version': xml_dictionary['general_info']['creator']['@VERSION'],
        'k_points_weights': [eigenvalue['k_point']['@weight'] for eigenvalue in xml_dictionary['output']['band_structure']['ks_energies']],
        'number_of_k_points': xml_dictionary['output']['band_structure']['nks'],
        'q_real_space': xml_dictionary['output']['algorithmic_info']['real_space_q'],
        'number_of_electrons': xml_dictionary['output']['band_structure']['nelec'],
    }

    if 'boundary_conditions' in xml_dictionary['output'] and 'assume_isolated' in xml_dictionary['output']['boundary_conditions']:
        xml_data['assume_isolated'] = xml_dictionary['output']['boundary_conditions']['assume_isolated']

    if 'fermi_energy' in xml_dictionary['output']['band_structure']:
        xml_data['fermi_energy'] = xml_dictionary['output']['band_structure']['fermi_energy']

    # We should put the `non_periodic_cell_correction` string in
    structure_data = {
       'atomic_positions_units': 'Angstrom',
       'direct_lattice_vectors_units': 'Angstrom',
        # ??? 'atoms_if_pos_list': [[1, 1, 1], [1, 1, 1]],
        'number_of_atoms': xml_dictionary['output']['atomic_structure']['@nat'],
        'lattice_parameter': xml_dictionary['input']['atomic_structure']['@alat'],
        'reciprocal_lattice_vectors': [
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b1'],
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b2'],
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b3']
        ],
        'atoms': [[atom['@name'], atom['$']] for atom in xml_dictionary['output']['atomic_structure']['atomic_positions']['atom']],
        'cell': {
            'lattice_vectors': lattice_vectors,
            'volume': cell_volume(*lattice_vectors),
            'atoms': [[atom['@name'], atom['$']] for atom in xml_dictionary['output']['atomic_structure']['atomic_positions']['atom']],
        },
        'lattice_parameter_xml': xml_dictionary['input']['atomic_structure']['@alat'],
        'number_of_species': xml_dictionary['input']['atomic_species']['@ntyp'],
        'species': {
            'index': [i + 1 for i,specie in enumerate(xml_dictionary['output']['atomic_species']['species'])],
            'pseudo': [specie['pseudo_file'] for specie in xml_dictionary['output']['atomic_species']['species']],
            'mass': [specie['mass'] for specie in xml_dictionary['output']['atomic_species']['species']],
            'type': [specie['@name'] for specie in xml_dictionary['output']['atomic_species']['species']]
        },
    }

    return xml_data, structure_data, {}
