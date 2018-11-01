# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import xmlschema
from xmlschema.exceptions import URLError
from defusedxml import ElementTree
import numpy as np

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.constants import ry_to_ev, hartree_to_ev, bohr_to_ang, ry_si, bohr_si, e_bohr2_to_coulomb_m2
from .versions import get_schema_filepath, get_default_schema_filepath


class QEXMLParsingError(QEOutputParsingError):
    pass


def raise_parsing_error(message):
    raise QEXMLParsingError(message)


def parser_assert(condition, message, log_func=raise_parsing_error):
    if not condition:
        log_func(message)


def parser_assert_equal(val1, val2, message, log_func=raise_parsing_error):
    if not (val1 == val2):
        msg = "Violated assertion: {} == {}".format(val1, val2)
        if message:
            msg += " - "
            msg += message
        log_func(msg)


def cell_volume(a1, a2, a3):
    r"""Returns the volume of the primitive cell: :math:`|\vec a_1\cdot(\vec a_2\cross \vec a_3)|`"""
    a_mid_0 = a2[1] * a3[2] - a2[2] * a3[1]
    a_mid_1 = a2[2] * a3[0] - a2[0] * a3[2]
    a_mid_2 = a2[0] * a3[1] - a2[1] * a3[0]

    return abs(float(a1[0] * a_mid_0 + a1[1] * a_mid_1 + a1[2] * a_mid_2))


def parse_pw_xml_post_6_2(xml_file, parser_opts, logger):
    """
    NOTE: the syntax for accessing the decoded XML dictionary (xml_dictionary)
    is the following:
          - If tag ['key'] is "simple", xml_dictionary['key'] returns its content;
          - otherwise:
            - xml_dictionary['key']['$'] returns its content
            - xml_dictionary['key']['@attr'] returns its attribute 'attr'
            - xml_dictionary['key']['nested_key'] goes one level deeper.
    """
    include_deprecated_v2_keys = parser_opts.get('include_deprecated_v2_keys', False)

    try:
        xml = ElementTree.parse(xml_file)
    except IOError:
        raise QEXMLParsingError('Could not open or parse the XML file {}'.format(xml_file))

    # detect schema name+path from XML contents
    schema_filepath = get_schema_filepath(xml)

    try:
        xsd = xmlschema.XMLSchema(schema_filepath)
    except URLError:

        # If loading the XSD file specified in the XML file fails, we try the default
        schema_filepath_default = get_default_schema_filepath()

        try:
            xsd = xmlschema.XMLSchema(schema_filepath_default)
        except URLError:
            raise QEXMLParsingError('Could not open or parse the XSD files {} and {}'.format(schema_filepath, schema_filepath_default))
        else:
            schema_filepath = schema_filepath_default

    # validate XML document against the schema
    xml_dictionary, errors = xsd.to_dict(xml, validation='lax')
    if errors:
        logger.error("{} XML schema validation error(s) (document: {} - schema: {}):".format(len(errors), xml_file, schema_filepath))
        for err in errors:
            logger.error(str(err))

    lattice_vectors = [
        [x * bohr_to_ang for x in xml_dictionary['output']['atomic_structure']['cell']['a1']],
        [x * bohr_to_ang for x in xml_dictionary['output']['atomic_structure']['cell']['a2']],
        [x * bohr_to_ang for x in xml_dictionary['output']['atomic_structure']['cell']['a3']],
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

    #if 'bands' in xml_dictionary['input'] and 'occupations' in xml_dictionary['input']['bands']:
    # the above condition is always true, according to the new schema

    try:
        occupations = xml_dictionary['input']['bands']['occupations']['$']  # also present as ['output']['band_structure']['occupations_kind']
    except TypeError:  # "string indices must be integers" -- might have attribute 'nspin'
        occupations = xml_dictionary['input']['bands']['occupations']

    if include_deprecated_v2_keys:

        # NOTE: these flags are not super clear
        # TODO: remove them and just use 'occupations' as a string, with possible values "smearing","tetrahedra", ...

        if occupations == 'from_input':  # False for 'fixed'  # TODO: does this make sense?
            fixed_occupations = True
        else:
            fixed_occupations = False

        if 'tetrahedra' in occupations:
            tetrahedron_method = True
        else:
            tetrahedron_method = False

        if occupations == 'from_input':
            smearing_method = False
        else:
            smearing_method = True

        if smearing_method:
            if 'smearing' not in (list(xml_dictionary['output']['band_structure'].keys()) + list(xml_dictionary['input']['bands'].keys())):
                logger.error("occupations is '{}' but key 'smearing' is not present under input/bands "
                             "nor under output/band_structure".format(occupations))
            # TODO: this error is triggered if no smearing is specified in input. But this is a valid input, so we should't throw an error.
            # Should we ask QE to print something nonetheless?
            # Also happens if occupations='fixed'.
            # (Example: occupations is 'fixed' but key 'smearing' is not present under input/bandsnor output/band_structure. See calculation 4940, versus 4981.)

    # TODO/NOTE: Not including smearing type and width for now.
    # In the old XML format they are under OCCUPATIONS as SMEARING_TYPE and SMEARING_PARAMETER,
    # but watch out: the value in the old format is half of that in the new format
    # (the code divides it by e2=2.0, see PW/src/pw_restart.f90:446)
    '''
    if 'smearing' in xml_dictionary['output']['band_structure']:
        smearing_xml = xml_dictionary['output']['band_structure']['smearing']
    elif 'smearing' in xml_dictionary['input']['bands']:
        smearing_xml = xml_dictionary['input']['bands']['smearing']
    try:
        smearing_type    = smearing_xml['$']
        smearing_degauss = smearing_xml['@degauss']
    except NameError:
        pass
    '''

    # Here are some notes from the QE code for reference.
    # SMEARING_METHOD = lgauss
    #       lgauss,         &! if .TRUE.: use gaussian broadening
    #       ltetra,         &! if .TRUE.: use tetrahedra
    # SMEARING_TYPE = ngauss  (see Modules/qexml.f90:1530)
    #       ngauss              ! type of smearing technique
    # From dos.x input description:
    #   Type of gaussian broadening:
    #      =  0  Simple Gaussian (default)
    #      =  1  Methfessel-Paxton of order 1
    #      = -1  Marzari-Vanderbilt "cold smearing"
    #      =-99  Fermi-Dirac function

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
    nsym = xml_dictionary['output']['symmetries']['nsym']  # crystal symmetries
    nrot = xml_dictionary['output']['symmetries']['nrot']  # lattice symmetries

    for symmetry in xml_dictionary['output']['symmetries']['symmetry']:

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
            sym['t_rev'] = u'1' if symmetry['info']['@time_reversal'] else u'0'
        except KeyError:
            sym['t_rev'] = u'0'

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
            raise QEXMLParsingError("Unexpected type of symmetry: {}".format(symmetry_type))

    if (nsym != len(symmetries)) or (nrot != len(symmetries)+len(lattice_symmetries)):
        logger.warning("Inconsistent number of symmetries: nsym={}, nrot={}, len(symmetries)={}, len(lattice_symmetries)={}"
                       .format(nsym, nrot, len(symmetries), len(lattice_symmetries))
        )

    # Band structure
    num_k_points   = xml_dictionary['output']['band_structure']['nks']
    num_electrons  = xml_dictionary['output']['band_structure']['nelec']
    num_atomic_wfc = xml_dictionary['output']['band_structure']['num_of_atomic_wfc']
    num_bands      = xml_dictionary['output']['band_structure'].get('nbnd')
    num_bands_up   = xml_dictionary['output']['band_structure'].get('nbnd_up')
    num_bands_down = xml_dictionary['output']['band_structure'].get('nbnd_dw')
    if (num_bands_up is None) and (num_bands_down is None):
        spins = False   # we are either in the non-polarized or non-collinear case
    else:
        # TODO: is it always nbnd_up==nbnd_dw ?
        parser_assert((num_bands_up is not None) and (num_bands_down is not None),
            "Only one of 'nbnd_up' and 'nbnd_dw' was found")
        if num_bands is not None:
            parser_assert(num_bands == num_bands_up + num_bands_down,
                "Inconsistent number of bands: nbnd={}, nbnd_up={}, nbnd_down={}".format(num_bands, num_bands_up, num_bands_down))
        else:
            num_bands = num_bands_up + num_bands_down   # backwards compatibility;
        spins = True

    # k-points
    k_points = []
    k_points_weights = []
    ks_states = xml_dictionary['output']['band_structure']['ks_energies']
    # alat is technically an optional attribute according to the schema,
    # but I don't know what to do if it's missing. atomic_structure is mandatory.
    output_alat_bohr = xml_dictionary['output']['atomic_structure']['@alat']
    output_alat_angstrom = output_alat_bohr * bohr_to_ang
    for ks_state in ks_states:
        k_points.append([kp*2*np.pi/output_alat_angstrom for kp in ks_state['k_point']['$']])
        k_points_weights.append(ks_state['k_point']['@weight'])
    # bands
    if not spins:
        # Note: 'parse_with_retrieved' still expects a list of lists
        band_eigenvalues = [[]]
        band_occupations = [[]]
        for ks_state in ks_states:
            band_eigenvalues[0].append(ks_state['eigenvalues']['$'])
            band_occupations[0].append(ks_state['occupations']['$'])
    else:
        band_eigenvalues = [[],[]]
        band_occupations = [[],[]]
        for ks_state in ks_states:
            band_eigenvalues[0].append(ks_state['eigenvalues']['$'][0:num_bands_up])
            band_eigenvalues[1].append(ks_state['eigenvalues']['$'][num_bands_up:num_bands])
            band_occupations[0].append(ks_state['occupations']['$'][0:num_bands_up])
            band_occupations[1].append(ks_state['occupations']['$'][num_bands_up:num_bands])

    band_eigenvalues = np.array(band_eigenvalues) * hartree_to_ev
    band_occupations = np.array(band_occupations)

    if not spins:
        parser_assert_equal(band_eigenvalues.shape, (1,num_k_points,num_bands),
                            "Unexpected shape of band_eigenvalues")
        parser_assert_equal(band_occupations.shape, (1,num_k_points,num_bands),
                            "Unexpected shape of band_occupations")
    else:
        parser_assert_equal(band_eigenvalues.shape, (2,num_k_points,num_bands_up),
                            "Unexpected shape of band_eigenvalues")
        parser_assert_equal(band_occupations.shape, (2,num_k_points,num_bands_up),
                            "Unexpected shape of band_occupations")

    bands_dict = {
        'occupations': band_occupations,
        'bands': band_eigenvalues,
        'bands_units': 'eV',
    }

    xml_data = {
        #'pp_check_flag': True, # Currently not printed in the new format.
            #Signals whether the XML file is complete
            # and can be used for post-processing. Everything should be in the XML now, but in
            # any case, the new XML schema should mostly protect from incomplete files.
        'lkpoint_dir': False, # Currently not printed in the new format.
            # Signals whether kpt-data are written in sub-directories.
            # Was generally true in the old format, but now all the eigenvalues are
            # in the XML file, under output / band_structure, so this is False.
        'charge_density': u'./charge-density.dat', # A file name. Not printed in the new format.
            # The filename and path are considered fixed: <outdir>/<prefix>.save/charge-density.dat
            # TODO: change to .hdf5 if output format is HDF5 (issue #222)
        'xml_warnings': [],
        'rho_cutoff_units': 'eV',
        'wfc_cutoff_units': 'eV',
        'fermi_energy_units': 'eV',
        'k_points_units': '1 / angstrom',
        'symmetries_units': 'crystal',
        'constraint_mag': constraint_mag,
        'occupations': occupations,
        'magnetization_angle2': magnetization_angle2,
        'magnetization_angle1': magnetization_angle1,
        'starting_magnetization': starting_magnetization,
        'has_electric_field': has_electric_field,
        'has_dipole_correction': has_dipole_correction,
        'lda_plus_u_calculation': 'dftU' in xml_dictionary['output'],
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
        'do_not_use_time_reversal': xml_dictionary['input']['symmetry_flags']['noinv'],
        'spin_orbit_domag': xml_dictionary['output']['magnetization']['do_magnetization'],
        'fft_grid': list(xml_dictionary['output']['basis_set']['fft_grid'].values()),
        'lsda': lsda,
        'number_of_spin_components': nspin,
        'no_time_rev_operations': xml_dictionary['input']['symmetry_flags']['no_t_rev'],
        'inversion_symmetry': inversion_symmetry,  # the old tag was INVERSION_SYMMETRY,
                    #and was set to (from the code): "invsym    if true the system has inversion symmetry"
        'number_of_bravais_symmetries': nrot,  # lattice symmetries
        'number_of_symmetries': nsym,          # crystal symmetries
        'wfc_cutoff': xml_dictionary['input']['basis']['ecutwfc'] * hartree_to_ev,
        'rho_cutoff': xml_dictionary['output']['basis_set']['ecutrho'] * hartree_to_ev, # not always printed in input->basis
        'smooth_fft_grid': list(xml_dictionary['output']['basis_set']['fft_smooth'].values()),
        'dft_exchange_correlation': xml_dictionary['input']['dft']['functional'],  # TODO: also parse optional elements of 'dft' tag
            # WARNING: this is different between old XML and new XML
        'spin_orbit_calculation': spin_orbit_calculation,
        'number_of_atomic_wfc': num_atomic_wfc,
        'number_of_k_points': num_k_points,
        'number_of_electrons': num_electrons,
        'k_points': k_points,
        'k_points_weights': k_points_weights,
        'q_real_space': xml_dictionary['output']['algorithmic_info']['real_space_q'],
    }

    try:
        monkhorst_pack = xml_dictionary['input']['k_points_IBZ']['monkhorst_pack']
    except KeyError:
        pass  # not using Monkhorst pack
    else:
        xml_data['monkhorst_pack_grid'] = [monkhorst_pack[attr] for attr in ['@nk1','@nk2','@nk3']]
        xml_data['monkhorst_pack_offset'] = [monkhorst_pack[attr] for attr in ['@k1','@k2','@k3']]

    if not spins:
        xml_data['number_of_bands'] = num_bands
    else:
        # QE counts them twice if spin-collinear
        # if non-collinear-spin, num_bands is None
        if num_bands % 2:
            xml_data['number_of_bands'] = float(num_bands)/2
        else:
            xml_data['number_of_bands'] = num_bands/2

    for key,value in [('number_of_bands_up',num_bands_up), ('number_of_bands_down',num_bands_down)]:
        if value is not None:
            xml_data[key] = value

    if 'boundary_conditions' in xml_dictionary['output'] and 'assume_isolated' in xml_dictionary['output']['boundary_conditions']:
        xml_data['assume_isolated'] = xml_dictionary['output']['boundary_conditions']['assume_isolated']

    if 'fermi_energy' in xml_dictionary['output']['band_structure']:
        xml_data['fermi_energy'] = xml_dictionary['output']['band_structure']['fermi_energy'] * hartree_to_ev

    # This is not printed by QE 6.3, but will be re-added before the next version
    if 'real_space_beta' in xml_dictionary['output']['algorithmic_info']:
        xml_data['beta_real_space'] = xml_dictionary['output']['algorithmic_info']['real_space_beta']

    conv_info = {}
    conv_info_scf = {}
    conv_info_opt = {}
    # NOTE: n_scf_steps refers to the number of SCF steps in the *last* loop only.
    # To get the total number of SCF steps in the run you should sum up the individual steps.
    # TODO: should we parse 'steps' too? Are they already added in the output trajectory?
    for key in ['convergence_achieved','n_scf_steps','scf_error']:
        try:
            conv_info_scf[key] = xml_dictionary['output']['convergence_info']['scf_conv'][key]
        except KeyError:
            pass
    for key in ['convergence_achieved','n_opt_steps','grad_norm']:
        try:
            conv_info_opt[key] = xml_dictionary['output']['convergence_info']['opt_conv'][key]
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

    if include_deprecated_v2_keys:
        xml_data['fixed_occupations'] = fixed_occupations
        xml_data['tetrahedron_method'] = tetrahedron_method
        xml_data['smearing_method'] = smearing_method

    try:
        berry_phase = xml_dictionary['output']['electric_field']['BerryPhase']
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
        xml_data['total_phase']      = berry_phase['totalPhase']['$']
        xml_data['total_phase_units']      = '2pi'
        xml_data['ionic_phase']      = berry_phase['totalPhase']['@ionic']
        xml_data['ionic_phase_units']      = '2pi'
        xml_data['electronic_phase'] = berry_phase['totalPhase']['@electronic']
        xml_data['electronic_phase_units'] = '2pi'
        polarization = berry_phase['totalPolarization']['polarization']['$']
        polarization_units = berry_phase['totalPolarization']['polarization']['@Units']
        polarization_modulus = berry_phase['totalPolarization']['modulus']
        parser_assert(polarization_units in ['e/bohr^2','C/m^2'],
                      "Unsupported units '{}' of total polarization".format(polarization_units))
        if polarization_units == 'e/bohr^2':
            polarization *= e_bohr2_to_coulomb_m2
            polarization_modulus *= e_bohr2_to_coulomb_m2
        xml_data['polarization']        = polarization
        xml_data['polarization_module'] = polarization_modulus # should be called "modulus"
        xml_data['polarization_units']  = 'C / m^2'
        xml_data['polarization_direction'] = berry_phase['totalPolarization']['direction']
        # TODO: add conversion for (e/Omega).bohr (requires to know Omega, the volume of the cell)
        # TODO (maybe): Not parsed:
        # - individual ionic phases
        # - individual electronic phases and weights

    # TODO: We should put the `non_periodic_cell_correction` string in (?)
    atoms = [[atom['@name'], [coord*bohr_to_ang for coord in atom['$']]] for atom in xml_dictionary['output']['atomic_structure']['atomic_positions']['atom']]
    species = xml_dictionary['output']['atomic_species']['species']
    structure_data = {
       'atomic_positions_units': 'Angstrom',
       'direct_lattice_vectors_units': 'Angstrom',
        # ??? 'atoms_if_pos_list': [[1, 1, 1], [1, 1, 1]],
        'number_of_atoms': xml_dictionary['output']['atomic_structure']['@nat'],
        'lattice_parameter': output_alat_angstrom,
        'reciprocal_lattice_vectors': [
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b1'],
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b2'],
            xml_dictionary['output']['basis_set']['reciprocal_lattice']['b3']
        ],
        'atoms': atoms,
        'cell': {
            'lattice_vectors': lattice_vectors,
            'volume': cell_volume(*lattice_vectors),
            'atoms': atoms,
        },
        'lattice_parameter_xml': output_alat_bohr,
        'number_of_species': xml_dictionary['output']['atomic_species']['@ntyp'],
        'species': {
            'index': [i + 1 for i,specie in enumerate(species)],
            'pseudo': [specie['pseudo_file'] for specie in species],
            'mass': [specie['mass'] for specie in species],
            'type': [specie['@name'] for specie in species]
        },
    }

    return xml_data, structure_data, bands_dict
