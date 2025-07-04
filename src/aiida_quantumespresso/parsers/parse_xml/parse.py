# -*- coding: utf-8 -*-
from urllib.error import URLError

import numpy as np
from packaging.version import Version
from qe_tools import CONSTANTS
from xmlschema import XMLSchema

from aiida_quantumespresso.utils.mapping import get_logging_container

from .exceptions import XMLParseError
from .versions import get_default_schema_filepath, get_schema_filepath


def raise_parsing_error(message):
    raise XMLParseError(message)


def parser_assert(condition, message, log_func=raise_parsing_error):
    if not condition:
        log_func(message)


def parser_assert_equal(val1, val2, message, log_func=raise_parsing_error):
    if not (val1 == val2):
        msg = f'Violated assertion: {val1} == {val2}'
        if message:
            msg += ' - '
            msg += message
        log_func(msg)


def cell_volume(a1, a2, a3):
    r"""Returns the volume of the primitive cell: :math:`|\vec a_1\cdot(\vec a_2\cross \vec a_3)|`"""
    a_mid_0 = a2[1] * a3[2] - a2[2] * a3[1]
    a_mid_1 = a2[2] * a3[0] - a2[0] * a3[2]
    a_mid_2 = a2[0] * a3[1] - a2[1] * a3[0]

    return abs(float(a1[0] * a_mid_0 + a1[1] * a_mid_1 + a1[2] * a_mid_2))


def parse_xml_post_6_2(xml):
    """Parse the content of XML output file written by `pw.x` and `cp.x` with the new schema-based XML format.

    :param xml: parsed XML
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    e_bohr2_to_coulomb_m2 = 57.214766  # e/a0^2 to C/m^2 (electric polarization) from Wolfram Alpha

    logs = get_logging_container()

    schema_filepath = get_schema_filepath(xml)

    try:
        xsd = XMLSchema(schema_filepath)
    except URLError:

        # If loading the XSD file specified in the XML file fails, we try the default
        schema_filepath_default = get_default_schema_filepath()

        try:
            xsd = XMLSchema(schema_filepath_default)
        except URLError:
            raise XMLParseError(
                f'Could not open or parse the XSD files {schema_filepath} and {schema_filepath_default}'
            )
        else:
            schema_filepath = schema_filepath_default

    # Validate XML document against the schema
    # Returned dictionary has a structure where, if tag ['key'] is "simple", xml_dictionary['key'] returns its content.
    # Otherwise, the following keys are available:
    #
    #  xml_dictionary['key']['$'] returns its content
    #  xml_dictionary['key']['@attr'] returns its attribute 'attr'
    #  xml_dictionary['key']['nested_key'] goes one level deeper.

    # Fix a bug of QE 6.8: the output XML is not consistent with schema, see
    # https://github.com/aiidateam/aiida-quantumespresso/pull/717
    xml_creator = xml.find('./general_info/creator')
    if xml_creator is not None and 'VERSION' in xml_creator.attrib:
        creator_version = xml_creator.attrib['VERSION']
        if creator_version == '6.8':
            root = xml.getroot()
            timing_info = root.find('./timing_info')
            partial_pwscf = timing_info.find("partial[@label='PWSCF'][@calls='0']")
            try:
                timing_info.remove(partial_pwscf)
            except (TypeError, ValueError):
                pass

    xml_dictionary, errors = xsd.to_dict(xml, validation='lax')
    if errors:
        logs.error.append(f'{len(errors)} XML schema validation error(s) schema: {schema_filepath}:')
        for err in errors:
            logs.error.append(str(err))

    xml_version = Version(xml_dictionary['general_info']['xml_format']['@VERSION'])
    inputs = xml_dictionary.get('input', {})
    outputs = xml_dictionary['output']

    lattice_vectors = [
        [x * CONSTANTS.bohr_to_ang for x in outputs['atomic_structure']['cell']['a1']],
        [x * CONSTANTS.bohr_to_ang for x in outputs['atomic_structure']['cell']['a2']],
        [x * CONSTANTS.bohr_to_ang for x in outputs['atomic_structure']['cell']['a3']],
    ]

    has_electric_field = inputs.get('electric_field', {}).get('electric_potential', None) == 'sawtooth_potential'
    has_dipole_correction = inputs.get('electric_field', {}).get('dipole_correction', False)

    if 'occupations' in inputs.get('bands', {}):
        try:
            occupations = inputs['bands']['occupations']['$']  # yapf: disable
        except TypeError:  # "string indices must be integers" -- might have attribute 'nspin'
            occupations = inputs['bands']['occupations']
    else:
        occupations = None

    starting_magnetization = []
    magnetization_angle1 = []
    magnetization_angle2 = []

    for specie in outputs['atomic_species']['species']:
        starting_magnetization.append(specie.get('starting_magnetization', 0.0))
        magnetization_angle1.append(specie.get('magnetization_angle1', 0.0))
        magnetization_angle2.append(specie.get('magnetization_angle2', 0.0))

    constraint_mag = 0
    spin_constraints = inputs.get('spin_constraints', {}).get('spin_constraints', None)
    if spin_constraints == 'atomic':
        constraint_mag = 1
    elif spin_constraints == 'atomic direction':
        constraint_mag = 2
    elif spin_constraints == 'total':
        constraint_mag = 3
    elif spin_constraints == 'total direction':
        constraint_mag = 6

    lsda = inputs.get('spin', {}).get('lsda', False)
    spin_orbit_calculation = inputs.get('spin', {}).get('spinorbit', False)
    non_colinear_calculation = outputs.get('magnetization', {}).get('noncolin', False)
    do_magnetization = outputs.get('magnetization', {}).get('do_magnetization', False)

    # Time reversal symmetry of the system
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

    symmetries = []
    lattice_symmetries = []  # note: will only contain lattice symmetries that are NOT crystal symmetries
    inversion_symmetry = False

    # See also PW/src/setup.f90
    nsym = outputs.get('symmetries', {}).get('nsym', None)  # crystal symmetries
    nrot = outputs.get('symmetries', {}).get('nrot', None)  # lattice symmetries

    for symmetry in outputs.get('symmetries', {}).get('symmetry', []):

        # There are two types of symmetries, lattice and crystal. The pure inversion (-I) is always a lattice symmetry,
        # so we don't care. But if the pure inversion is also a crystal symmetry, then then the system as a whole
        # has (by definition) inversion symmetry; so we set the global property inversion_symmetry = True.
        symmetry_type = symmetry['info']['$']
        symmetry_name = symmetry['info']['@name']
        if symmetry_type == 'crystal_symmetry' and symmetry_name.lower() == 'inversion':
            inversion_symmetry = True

        sym = {
            'rotation': [
                symmetry['rotation']['$'][0:3],
                symmetry['rotation']['$'][3:6],
                symmetry['rotation']['$'][6:9],
            ],
            'name': symmetry_name,
        }

        try:
            sym['t_rev'] = '1' if symmetry['info']['@time_reversal'] else '0'
        except KeyError:
            sym['t_rev'] = '0'

        try:
            sym['equivalent_atoms'] = symmetry['equivalent_atoms']['$']
        except KeyError:
            pass

        try:
            sym['fractional_translation'] = symmetry['fractional_translation']
        except KeyError:
            pass

        if symmetry_type == 'crystal_symmetry':
            symmetries.append(sym)
        elif symmetry_type == 'lattice_symmetry':
            lattice_symmetries.append(sym)
        else:
            raise XMLParseError(f'Unexpected type of symmetry: {symmetry_type}')

    if (nsym != len(symmetries)) or (nrot != len(symmetries) + len(lattice_symmetries)):
        logs.warning.append(
            'Inconsistent number of symmetries: nsym={}, nrot={}, len(symmetries)={}, len(lattice_symmetries)={}'.
            format(nsym, nrot, len(symmetries), len(lattice_symmetries))
        )

    xml_data = {
        #'pp_check_flag': True, # Currently not printed in the new format.
        # Signals whether the XML file is complete
        # and can be used for post-processing. Everything should be in the XML now, but in
        # any case, the new XML schema should mostly protect from incomplete files.
        'lkpoint_dir': False,  # Currently not printed in the new format.
        # Signals whether kpt-data are written in sub-directories.
        # Was generally true in the old format, but now all the eigenvalues are
        # in the XML file, under output / band_structure, so this is False.
        'charge_density': './charge-density.dat',  # A file name. Not printed in the new format.
        # The filename and path are considered fixed: <outdir>/<prefix>.save/charge-density.dat
        # TODO: change to .hdf5 if output format is HDF5 (issue #222)
        'rho_cutoff_units': 'eV',
        'wfc_cutoff_units': 'eV',
        'fermi_energy_units': 'eV',
        'k_points_units': '1 / angstrom',
        'symmetries_units': 'crystal',
        'constraint_mag': constraint_mag,
        'magnetization_angle2': magnetization_angle2,
        'magnetization_angle1': magnetization_angle1,
        'starting_magnetization': starting_magnetization,
        'has_electric_field': has_electric_field,
        'has_dipole_correction': has_dipole_correction,
        'lda_plus_u_calculation': 'dftU' in outputs,
        'format_name': xml_dictionary['general_info']['xml_format']['@NAME'],
        'format_version': xml_dictionary['general_info']['xml_format']['@VERSION'],
        # TODO: check that format version: a) matches the XSD schema version; b) is updated as well
        #       See line 43 in Modules/qexsd.f90
        'creator_name': xml_dictionary['general_info']['creator']['@NAME'].lower(),
        'creator_version': xml_dictionary['general_info']['creator']['@VERSION'],
        'non_colinear_calculation': non_colinear_calculation,
        'do_magnetization': do_magnetization,
        'time_reversal_flag': time_reversal,
        'symmetries': symmetries,
        'lattice_symmetries': lattice_symmetries,
        'do_not_use_time_reversal': inputs.get('symmetry_flags', {}).get('noinv', None),
        'spin_orbit_domag': do_magnetization,
        'fft_grid': [value for _, value in sorted(outputs['basis_set']['fft_grid'].items())],
        'lsda': lsda,
        'number_of_spin_components': nspin,
        'no_time_rev_operations': inputs.get('symmetry_flags', {}).get('no_t_rev', None),
        'inversion_symmetry':
        inversion_symmetry,  # the old tag was INVERSION_SYMMETRY and was set to (from the code): "invsym    if true the system has inversion symmetry"
        'number_of_bravais_symmetries': nrot,  # lattice symmetries
        'number_of_symmetries': nsym,  # crystal symmetries
        'wfc_cutoff': inputs.get('basis', {}).get('ecutwfc', -1.0) * CONSTANTS.hartree_to_ev,
        'rho_cutoff': outputs['basis_set']['ecutrho'] * CONSTANTS.hartree_to_ev,  # not always printed in input->basis
        'smooth_fft_grid': [value for _, value in sorted(outputs['basis_set']['fft_smooth'].items())],
        'dft_exchange_correlation': inputs.get('dft', {}).get('functional',
                                                              None),  # TODO: also parse optional elements of 'dft' tag
        # WARNING: this is different between old XML and new XML
        'spin_orbit_calculation': spin_orbit_calculation,
        'q_real_space': outputs['algorithmic_info']['real_space_q'],
    }

    # alat is technically an optional attribute according to the schema,
    # but I don't know what to do if it's missing. atomic_structure is mandatory.
    output_alat_bohr = outputs['atomic_structure']['@alat']
    output_alat_angstrom = output_alat_bohr * CONSTANTS.bohr_to_ang

    # Band structure
    if 'band_structure' in outputs:
        band_structure = outputs['band_structure']

        smearing_xml = None

        if 'smearing' in outputs['band_structure']:
            smearing_xml = outputs['band_structure']['smearing']
        elif 'smearing' in inputs:
            smearing_xml = inputs['smearing']

        if smearing_xml:
            degauss = smearing_xml['@degauss']

            # Versions below 19.03.04 (Quantum ESPRESSO<=6.4.1) incorrectly print degauss in Ry instead of Hartree
            if xml_version < Version('19.03.04'):
                degauss *= CONSTANTS.ry_to_ev
            else:
                degauss *= CONSTANTS.hartree_to_ev

            xml_data['degauss'] = degauss
            xml_data['smearing_type'] = smearing_xml['$']

        num_k_points = band_structure['nks']
        num_electrons = band_structure['nelec']

        # In schema v240411 (QE v7.3.1), the `number_of_atomic_wfc` is moved to the `atomic_structure` tag as an attribute
        num_atomic_wfc = band_structure.get('num_of_atomic_wfc', None) or outputs['atomic_structure']['@num_of_atomic_wfc']
        num_bands = band_structure.get('nbnd', None)
        num_bands_up = band_structure.get('nbnd_up', None)
        num_bands_down = band_structure.get('nbnd_dw', None)

        if num_bands is None and num_bands_up is None and num_bands_down is None:
            raise XMLParseError('None of `nbnd`, `nbnd_up` or `nbdn_dw` could be parsed.')

        # If both channels are `None` we are dealing with a non spin-polarized or non-collinear calculation
        elif num_bands_up is None and num_bands_down is None:
            spins = False

        # If only one of the channels is `None` we raise, because that is an inconsistent result
        elif num_bands_up is None or num_bands_down is None:
            raise XMLParseError('Only one of `nbnd_up` and `nbnd_dw` could be parsed')

        # Here it is a spin-polarized calculation, where for pw.x the number of bands in each channel should be identical.
        else:
            spins = True
            if num_bands_up != num_bands_down:
                raise XMLParseError(f'different number of bands for spin channels: {num_bands_up} and {num_bands_down}')

            if num_bands is not None and num_bands != num_bands_up + num_bands_down:
                raise XMLParseError(
                    'Inconsistent number of bands: nbnd={}, nbnd_up={}, nbnd_down={}'.format(
                        num_bands, num_bands_up, num_bands_down
                    )
                )

            if num_bands is None:
                num_bands = num_bands_up + num_bands_down  # backwards compatibility;

        k_points = []
        k_points_weights = []
        ks_states = band_structure['ks_energies']
        for ks_state in ks_states:
            k_points.append([kp * 2 * np.pi / output_alat_angstrom for kp in ks_state['k_point']['$']])
            k_points_weights.append(ks_state['k_point']['@weight'])

        if not spins:
            band_eigenvalues = [[]]
            band_occupations = [[]]
            for ks_state in ks_states:
                band_eigenvalues[0].append(ks_state['eigenvalues']['$'])
                band_occupations[0].append(ks_state['occupations']['$'])
        else:
            band_eigenvalues = [[], []]
            band_occupations = [[], []]
            for ks_state in ks_states:
                band_eigenvalues[0].append(ks_state['eigenvalues']['$'][0:num_bands_up])
                band_eigenvalues[1].append(ks_state['eigenvalues']['$'][num_bands_up:num_bands])
                band_occupations[0].append(ks_state['occupations']['$'][0:num_bands_up])
                band_occupations[1].append(ks_state['occupations']['$'][num_bands_up:num_bands])

        band_eigenvalues = np.array(band_eigenvalues) * CONSTANTS.hartree_to_ev
        band_occupations = np.array(band_occupations)

        if not spins:
            parser_assert_equal(
                band_eigenvalues.shape, (1, num_k_points, num_bands), 'Unexpected shape of band_eigenvalues'
            )
            parser_assert_equal(
                band_occupations.shape, (1, num_k_points, num_bands), 'Unexpected shape of band_occupations'
            )
        else:
            parser_assert_equal(
                band_eigenvalues.shape, (2, num_k_points, num_bands_up), 'Unexpected shape of band_eigenvalues'
            )
            parser_assert_equal(
                band_occupations.shape, (2, num_k_points, num_bands_up), 'Unexpected shape of band_occupations'
            )

        if not spins:
            xml_data['number_of_bands'] = num_bands
        else:
            # For collinear spin-polarized calculations `spins=True` and `num_bands` is sum of both channels. To get the
            # actual number of bands, we divide by two using integer division
            xml_data['number_of_bands'] = num_bands // 2

        for key, value in [('number_of_bands_up', num_bands_up), ('number_of_bands_down', num_bands_down)]:
            if value is not None:
                xml_data[key] = value

        if 'fermi_energy' in band_structure:
            xml_data['fermi_energy'] = band_structure['fermi_energy'] * CONSTANTS.hartree_to_ev

        if 'two_fermi_energies' in band_structure:
            xml_data['fermi_energy_up'], xml_data['fermi_energy_down'] = [
                energy * CONSTANTS.hartree_to_ev for energy in band_structure['two_fermi_energies']
            ]

        bands_dict = {
            'occupations': band_occupations,
            'bands': band_eigenvalues,
            'bands_units': 'eV',
        }

        xml_data['number_of_atomic_wfc'] = num_atomic_wfc
        xml_data['number_of_k_points'] = num_k_points
        xml_data['number_of_electrons'] = num_electrons
        xml_data['k_points'] = k_points
        xml_data['k_points_weights'] = k_points_weights
        xml_data['bands'] = bands_dict

    try:
        monkhorst_pack = inputs['k_points_IBZ']['monkhorst_pack']
    except KeyError:
        pass  # not using Monkhorst pack
    else:
        xml_data['monkhorst_pack_grid'] = [monkhorst_pack[attr] for attr in ['@nk1', '@nk2', '@nk3']]
        xml_data['monkhorst_pack_offset'] = [monkhorst_pack[attr] for attr in ['@k1', '@k2', '@k3']]

    if occupations is not None:
        xml_data['occupations'] = occupations

    if 'boundary_conditions' in outputs and 'assume_isolated' in outputs['boundary_conditions']:
        xml_data['assume_isolated'] = outputs['boundary_conditions']['assume_isolated']

    # This is not printed by QE 6.3, but will be re-added before the next version
    if 'real_space_beta' in outputs['algorithmic_info']:
        xml_data['beta_real_space'] = outputs['algorithmic_info']['real_space_beta']

    conv_info = {}
    conv_info_scf = {}
    conv_info_opt = {}
    # NOTE: n_scf_steps refers to the number of SCF steps in the *last* loop only.
    # To get the total number of SCF steps in the run you should sum up the individual steps.
    # TODO: should we parse 'steps' too? Are they already added in the output trajectory?
    for key in ['convergence_achieved', 'n_scf_steps', 'scf_error']:
        try:
            conv_info_scf[key] = outputs['convergence_info']['scf_conv'][key]
        except KeyError:
            pass
    for key in ['convergence_achieved', 'n_opt_steps', 'grad_norm']:
        try:
            conv_info_opt[key] = outputs['convergence_info']['opt_conv'][key]
        except KeyError:
            pass
    if conv_info_scf:
        conv_info['scf_conv'] = conv_info_scf
    if conv_info_opt:
        conv_info['opt_conv'] = conv_info_opt
    if conv_info:
        xml_data['convergence_info'] = conv_info

    if 'status' in xml_dictionary:
        xml_data['exit_status'] = xml_dictionary['status']
        # 0 = convergence reached;
        # -1 = SCF convergence failed;
        # 3 = ionic convergence failed
        # These might be changed in the future. Also see PW/src/run_pwscf.f90

    try:
        berry_phase = outputs['electric_field']['BerryPhase']
    except KeyError:
        pass
    else:
        # This is what I would like to do, but it's not retro-compatible
        # xml_data['berry_phase'] = {}
        # xml_data['berry_phase']['total_phase']         = berry_phase['totalPhase']['$']
        # xml_data['berry_phase']['total_phase_modulus'] = berry_phase['totalPhase']['@modulus']
        # xml_data['berry_phase']['total_ionic_phase']      = berry_phase['totalPhase']['@ionic']
        # xml_data['berry_phase']['total_electronic_phase'] = berry_phase['totalPhase']['@electronic']
        # xml_data['berry_phase']['total_polarization']           = berry_phase['totalPolarization']['polarization']['$']
        # xml_data['berry_phase']['total_polarization_modulus']   = berry_phase['totalPolarization']['modulus']
        # xml_data['berry_phase']['total_polarization_units']     = berry_phase['totalPolarization']['polarization']['@Units']
        # xml_data['berry_phase']['total_polarization_direction'] = berry_phase['totalPolarization']['direction']
        # parser_assert_equal(xml_data['berry_phase']['total_phase_modulus'].lower(), '(mod 2)',
        #                    "Unexpected modulus for total phase")
        # parser_assert_equal(xml_data['berry_phase']['total_polarization_units'].lower(), 'e/bohr^2',
        #                    "Unsupported units for total polarization")
        # Retro-compatible keys:
        polarization = berry_phase['totalPolarization']['polarization']['$']
        polarization_units = berry_phase['totalPolarization']['polarization']['@Units']
        polarization_modulus = berry_phase['totalPolarization']['modulus']
        parser_assert(
            polarization_units in ['e/bohr^2', 'C/m^2'],
            f"Unsupported units '{polarization_units}' of total polarization"
        )
        if polarization_units == 'e/bohr^2':
            polarization *= e_bohr2_to_coulomb_m2
            polarization_modulus *= e_bohr2_to_coulomb_m2

        xml_data['total_phase'] = berry_phase['totalPhase']['$']
        xml_data['total_phase_units'] = '2pi'
        xml_data['ionic_phase'] = berry_phase['totalPhase']['@ionic']
        xml_data['ionic_phase_units'] = '2pi'
        xml_data['electronic_phase'] = berry_phase['totalPhase']['@electronic']
        xml_data['electronic_phase_units'] = '2pi'
        xml_data['polarization'] = polarization
        xml_data['polarization_module'] = polarization_modulus  # should be called "modulus"
        xml_data['polarization_units'] = 'C / m^2'
        xml_data['polarization_direction'] = berry_phase['totalPolarization']['direction']
        # TODO: add conversion for (e/Omega).bohr (requires to know Omega, the volume of the cell)
        # TODO (maybe): Not parsed:
        # - individual ionic phases
        # - individual electronic phases and weights

    # TODO: We should put the `non_periodic_cell_correction` string in (?)
    atoms = [[atom['@name'], [coord * CONSTANTS.bohr_to_ang
                              for coord in atom['$']]]
             for atom in outputs['atomic_structure']['atomic_positions']['atom']]
    species = outputs['atomic_species']['species']
    structure_data = {
        'atomic_positions_units':
        'Angstrom',
        'direct_lattice_vectors_units':
        'Angstrom',
        # ??? 'atoms_if_pos_list': [[1, 1, 1], [1, 1, 1]],
        'number_of_atoms':
        outputs['atomic_structure']['@nat'],
        'lattice_parameter':
        output_alat_angstrom,
        'reciprocal_lattice_vectors': [
            outputs['basis_set']['reciprocal_lattice']['b1'], outputs['basis_set']['reciprocal_lattice']['b2'],
            outputs['basis_set']['reciprocal_lattice']['b3']
        ],
        'atoms':
        atoms,
        'cell': {
            'lattice_vectors': lattice_vectors,
            'volume': cell_volume(*lattice_vectors),
            'atoms': atoms,
        },
        'lattice_parameter_xml':
        output_alat_bohr,
        'number_of_species':
        outputs['atomic_species']['@ntyp'],
        'species': {
            'index': [i + 1 for i, specie in enumerate(species)],
            'pseudo': [specie['pseudo_file'] for specie in species],
            'mass': [specie['mass'] for specie in species],
            'type': [specie['@name'] for specie in species]
        },
    }

    xml_data['structure'] = structure_data

    return xml_data, logs
