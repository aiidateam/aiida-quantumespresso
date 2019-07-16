# -*- coding: utf-8 -*-
"""
A collection of function that are used to parse the output of Quantum Espresso PW.
The function that needs to be called from outside is parse_raw_output().
The functions mostly work without aiida specific functionalities.
The parsing will try to convert whatever it can in some dictionary, which
by operative decision doesn't have much structure encoded, [the values are simple ]
"""
from __future__ import absolute_import
from __future__ import print_function

import re
from six.moves import range

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.constants import ry_to_ev, bohr_to_ang, ry_si, bohr_si
from aiida_quantumespresso.utils.mapping import get_logging_container

lattice_tolerance = 1.e-5
units_suffix = '_units'
default_charge_units = 'e'
default_dipole_units = 'Debye'
default_energy_units = 'eV'
default_force_units = 'ev / angstrom'
default_k_points_units = '1 / angstrom'
default_length_units = 'Angstrom'
default_magnetization_units = 'Bohrmag / cell'
default_polarization_units = 'C / m^2'
default_stress_units = 'GPascal'


def detect_important_message(logs, line):

    message_map = {
        'error': {
            'Maximum CPU time exceeded': 'ERROR_OUT_OF_WALLTIME',
            'convergence NOT achieved after': 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            'history already reset at previous step: stopping': 'ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE',
            'problems computing cholesky': 'ERROR_DIAGONALIZATION_CHOLESKY_DECOMPOSITION',
            'charge is wrong': 'ERROR_CHARGE_IS_WRONG',
            'not orthogonal operation': 'ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION',
        },
        'warning': {
            'Warning:': None,
            'DEPRECATED:': None,
        }
    }

    # Match any known error and warning messages
    for marker, message in message_map['error'].items():
        if marker in line:
            if message is None:
                message = line
            logs.error.append(message)

    for marker, message in message_map['warning'].items():
        if marker in line:
            if message is None:
                message = line
            logs.warning.append(message)


def parse_stdout(stdout, input_parameters, parser_options=None, parsed_xml=None):
    """Parses the stdout content of a Quantum ESPRESSO `pw.x` calculation.

    :param stdout: the stdout content as a string
    :param input_parameters: dictionary with the input parameters
    :param parser_options: the parser options from the settings input parameter node
    :param parsed_xml: dictionary with data parsed from the XML output file
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    if parser_options is None:
        parser_options = {}

    if parsed_xml is None:
        parsed_xml = {}

    # Separate the input string into separate lines
    data_lines = stdout.split('\n')

    logs = get_logging_container()

    parsed_data = {}
    vdw_correction = False
    bands_data = parsed_xml.pop('bands', {})
    structure_data = parsed_xml.pop('structure', {})
    trajectory_data = {}

    maximum_ionic_steps = None
    marker_bfgs_converged = False

    # First check whether the `JOB DONE` message was written, otherwise the job was interrupted
    for line in data_lines:
        if 'JOB DONE' in line:
            break
    else:
        logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

    # Determine whether the input switched on an electric field
    lelfield = input_parameters.get('CONTROL', {}).get('lelfield', False)

    # Find some useful quantities.
    if not parsed_xml.get('number_of_bands', None) and not structure_data:
        try:
            for line in stdout.split('\n'):
                if 'lattice parameter (alat)' in line:
                    alat = float(line.split('=')[1].split('a.u')[0])
                elif 'number of atoms/cell' in line:
                    nat = int(line.split('=')[1])
                elif 'number of atomic types' in line:
                    ntyp = int(line.split('=')[1])
                elif 'unit-cell volume' in line:
                    if '(a.u.)^3' in line:
                        volume = float(line.split('=')[1].split('(a.u.)^3')[0])
                    else:
                        # occurs in v5.3.0
                        volume = float(line.split('=')[1].split('a.u.^3')[0])
                elif 'number of Kohn-Sham states' in line:
                    nbnd = int(line.split('=')[1])
                elif "number of k points" in line:
                    nk = int(line.split('=')[1].split()[0])
                    if input_parameters.get('SYSTEM', {}).get('nspin', 1) > 1:
                        # QE counts twice each k-point in spin-polarized calculations
                        nk /= 2
                elif "Dense  grid" in line:
                    FFT_grid = [int(g) for g in
                                line.split('(')[1].split(')')[0].split(',')]
                elif "Smooth grid" in line:
                    smooth_FFT_grid = [int(g) for g in
                                       line.split('(')[1].split(')')[0].split(',')]
                    break
            alat *= bohr_to_ang
            volume *= bohr_to_ang**3
            parsed_data['lattice_parameter_initial'] = alat
            parsed_data['number_of_bands'] = nbnd
            try:
                parsed_data['number_of_k_points'] = nk
                parsed_data['fft_grid'] = FFT_grid
                parsed_data['smooth_fft_grid'] = smooth_FFT_grid
            except NameError:  # these are not crucial, so parsing does not fail if they are not found
                pass
        except NameError:  # nat or other variables where not found, and thus not initialized

            # Try to get some error messages
            lines = stdout.split('\n')

            for line_number, line in enumerate(lines):
                # Compare the line to the known set of error and warning messages and add them to the log container
                detect_important_message(logs, line)

            if len(logs.error) or len(logs.warning) > 0:
                parsed_data['trajectory'] = trajectory_data
                return parsed_data, logs

            # did not find any error message -> raise an Error and do not return anything
            raise QEOutputParsingError('Parser cannot load basic info.')
    else:
        nat = structure_data['number_of_atoms']
        ntyp = structure_data['number_of_species']
        nbnd = parsed_xml['number_of_bands']
        alat = structure_data['lattice_parameter_xml']
        volume = structure_data['cell']['volume']
    # NOTE: lattice_parameter_xml is the lattice parameter of the xml file
    # in the units used by the code. lattice_parameter instead in angstroms.

    # Save these two quantities in the parsed_data, because they will be
    # useful for queries (maybe), and structure_data will not be stored as a Dict
    parsed_data['number_of_atoms'] = nat
    parsed_data['number_of_species'] = ntyp
    parsed_data['volume'] = volume

    c_bands_error = False

    # now grep quantities that can be considered isolated informations.
    for count, line in enumerate(data_lines):

        # Compare the line to the known set of error and warning messages and add them to the log container
        detect_important_message(logs, line)

        # to be used for later
        if 'Carrying out vdW-DF run using the following parameters:' in line:
            vdw_correction = True

        elif 'Cartesian axes' in line:
            # this is the part when initial positions and chemical
            # symbols are printed (they do not change during a run)
            i = count + 1
            while i < count + 10 and not('site n.' in data_lines[i] and
                                      'atom' in data_lines[i]):
                i += 1
            if 'site n.' in data_lines[i] and 'atom' in data_lines[i]:
                trajectory_data['atomic_species_name'] = [data_lines[i + 1 + j].split()[1] for j in range(nat)]

        # parse the initialization time (take only first occurence)
        elif ('init_wall_time_seconds' not in parsed_data and "total cpu time spent up to now is" in line):
            init_time = float(line.split("total cpu time spent up to now is"
                                         )[1].split('secs')[0])
            parsed_data['init_wall_time_seconds'] = init_time

        # parse the global file, for informations that are written only once
        elif 'PWSCF' in line and 'WALL' in line:
            try:
                time = line.split('CPU')[1].split('WALL')[0]
                parsed_data['wall_time'] = time
            except Exception:
                logs.warning.append('Error while parsing wall time.')
            try:
                parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
            except ValueError:
                raise QEOutputParsingError("Unable to convert wall_time in seconds.")

        # for later control on relaxation-dynamics convergence
        elif 'nstep' in line and '=' in line:
            maximum_ionic_steps = int(line.split()[2])

        elif 'bfgs converged in' in line:
            marker_bfgs_converged = True

        elif 'number of bfgs steps' in line:
            try:
                parsed_data['number_ionic_steps'] += 1
            except KeyError:
                parsed_data['number_ionic_steps'] = 1

        elif 'A final scf calculation at the relaxed structure' in line:
            parsed_data['final_scf'] = True

        elif 'point group' in line:
            if 'k-point group' not in line:
                try:
                    # Split line in components delimited by either space(s) or
                    # parenthesis and filter out empty strings
                    line_elems = [_f for _f in re.split(r' +|\(|\)', line) if _f]

                    pg_international = line_elems[-1]
                    pg_schoenflies = line_elems[-2]

                    parsed_data['pointgroup_international'] = pg_international
                    parsed_data['pointgroup_schoenflies'] = pg_schoenflies

                except Exception:
                    warning = "Problem parsing point group, I found: {}".format(line.strip())
                    logs.warning.append(warning)

        # special parsing of c_bands error
        elif 'c_bands' in line and 'eigenvalues not converged' in line:
            c_bands_error = True

        elif "iteration #" in line:
            if 'Calculation restarted' not in line and 'Calculation stopped' not in line:
                try:
                    parsed_data['total_number_of_scf_iterations'] += 1
                except KeyError:
                    parsed_data['total_number_of_scf_iterations'] = 1

            if c_bands_error:
                # if there is another iteration, c_bands is not necessarily a problem
                # I put a warning only if c_bands error appears in the last iteration
                c_bands_error = False

    if c_bands_error:
        logs.warning.append("c_bands: at least 1 eigenvalues not converged")

    # I split the output text in the atomic SCF calculations.
    # the initial part should be things already contained in the xml.
    # (cell, initial positions, kpoints, ...) and I skip them.
    # In case, parse for them before this point.
    # Put everything in a trajectory_data dictionary
    relax_steps = stdout.split('Self-consistent Calculation')[1:]
    relax_steps = [i.split('\n') for i in relax_steps]

    # now I create a bunch of arrays for every step.
    for data_step in relax_steps:

        trajectory_frame = {}

        for count, line in enumerate(data_step):

            if 'CELL_PARAMETERS' in line:
                try:
                    a1 = [float(s) for s in data_step[count + 1].split()]
                    a2 = [float(s) for s in data_step[count + 2].split()]
                    a3 = [float(s) for s in data_step[count + 3].split()]
                    # try except indexerror for not enough lines
                    lattice = line.split('(')[1].split(')')[0].split('=')
                    if lattice[0].lower() not in ['alat', 'bohr', 'angstrom']:
                        raise QEOutputParsingError('Error while parsing cell_parameters: ' +
                                                   'unsupported units {}'.format(lattice[0]))

                    if 'alat' in lattice[0].lower():
                        a1 = [alat * bohr_to_ang * float(s) for s in a1]
                        a2 = [alat * bohr_to_ang * float(s) for s in a2]
                        a3 = [alat * bohr_to_ang * float(s) for s in a3]
                        lattice_parameter_b = float(lattice[1])
                        if abs(lattice_parameter_b - alat) > lattice_tolerance:
                            raise QEOutputParsingError("Lattice parameters mismatch! " +
                                                       "{} vs {}".format(lattice_parameter_b, alat))
                    elif 'bohr' in lattice[0].lower():
                        lattice_parameter_b *= bohr_to_ang
                        a1 = [bohr_to_ang * float(s) for s in a1]
                        a2 = [bohr_to_ang * float(s) for s in a2]
                        a3 = [bohr_to_ang * float(s) for s in a3]
                    try:
                        trajectory_data['lattice_vectors_relax'].append([a1, a2, a3])
                    except KeyError:
                        trajectory_data['lattice_vectors_relax'] = [[a1, a2, a3]]

                except Exception:
                    logs.warning.append('Error while parsing relaxation cell parameters.')

            elif 'ATOMIC_POSITIONS' in line:
                try:
                    this_key = 'atomic_positions_relax'
                    # the inizialization of tau prevent parsed_data to be associated
                    # to the pointer of the previous iteration
                    metric = line.split('(')[1].split(')')[0]
                    if metric not in ['alat', 'bohr', 'angstrom']:
                        raise QEOutputParsingError('Error while parsing atomic_positions: units not supported.')
                    # TODO: check how to map the atoms in the original scheme
                    positions = []
                    for i in range(nat):
                        line2 = data_step[count + 1 + i].split()
                        tau = [float(s) for s in line2[1:4]]
                        if metric == 'alat':
                            tau = [alat * float(s) for s in tau]
                        elif metric == 'bohr':
                            tau = [bohr_to_ang * float(s) for s in tau]
                        positions.append(tau)
                    try:
                        trajectory_data[this_key].append(positions)
                    except KeyError:
                        trajectory_data[this_key] = [positions]
                except Exception:
                    logs.warning.append('Error while parsing relaxation atomic positions.')

            # NOTE: in the above, the chemical symbols are not those of AiiDA
            # since the AiiDA structure is different. So, I assume now that the
            # order of atoms is the same of the input atomic structure.

            # Computed dipole correction in slab geometries.
            # save dipole in debye units, only at last iteration of scf cycle
            elif 'Computed dipole along edir' in line:
                j = count + 3
                line2 = data_step[j]
                try:
                    units = line2.split()[-1]
                    if default_dipole_units.lower() not in units.lower():  # only debye
                        raise QEOutputParsingError("Error parsing the dipole correction. Units {} are not supported.".format(units))
                    value = float(line2.split()[-2])
                except IndexError:  # on units
                    pass
                # save only the last dipole correction
                while 'Computed dipole along edir' not in line2:
                    j += 1
                    try:
                        line2 = data_step[j]
                    except IndexError:  # The dipole is also written at the beginning of a new bfgs iteration
                        break
                    if 'End of self-consistent calculation' in line2:
                        try:
                            trajectory_data['dipole'].append(value)
                        except KeyError:
                            trajectory_data['dipole'] = [value]
                        parsed_data['dipole' + units_suffix] = default_dipole_units
                        break

            elif 'convergence has been achieved in' in line or 'convergence NOT achieved after' in line:
                try:
                    scf_iterations = int(line.split("iterations")[0].split()[-1])
                    try:
                        trajectory_data['scf_iterations'].append(scf_iterations)
                    except KeyError:
                        trajectory_data['scf_iterations'] = [scf_iterations]
                except Exception:
                    logs.warning.append('Error while parsing scf iterations.')

            elif 'End of self-consistent calculation' in line:
                # parse energy threshold for diagonalization algorithm
                try:
                    j = 0
                    while True:
                        j -= 1
                        line2 = data_step[count + j]
                        if 'ethr' in line2:
                            value = float(line2.split('=')[1].split(',')[0])
                            break
                    try:
                        trajectory_data['energy_threshold'].append(value)
                    except KeyError:
                        trajectory_data['energy_threshold'] = [value]
                except Exception:
                    logs.warning.append('Error while parsing ethr.')

                # parse final magnetic moments, if present
                try:
                    j = 0
                    while True:
                        j -= 1
                        line2 = data_step[count + j]
                        if 'Magnetic moment per site' in line2:
                            break
                        if 'iteration' in line2:
                            raise QEOutputParsingError
                    mag_moments = []
                    charges = []
                    while True:
                        j += 1
                        line2 = data_step[count + j]
                        if 'atom:' in line2:
                            mag_moments.append(float(line2.split('magn:')[1].split()[0]))
                            charges.append(float(line2.split('charge:')[1].split()[0]))
                        if len(mag_moments) == nat:
                            break
                    try:
                        trajectory_data['atomic_magnetic_moments'].append(mag_moments)
                        trajectory_data['atomic_charges'].append(charges)
                    except KeyError:
                        trajectory_data['atomic_magnetic_moments'] = [mag_moments]
                        trajectory_data['atomic_charges'] = [charges]
                    parsed_data['atomic_magnetic_moments' + units_suffix] = default_magnetization_units
                    parsed_data['atomic_charges' + units_suffix] = default_charge_units
                except QEOutputParsingError:
                    pass

            # grep energy and possibly, magnetization
            elif '!' in line:
                try:
                    for key in ['energy', 'energy_accuracy']:
                        if key not in trajectory_data:
                            trajectory_data[key] = []

                    En = float(line.split('=')[1].split('Ry')[0]) * ry_to_ev
                    E_acc = float(data_step[count + 2].split('<')[1].split('Ry')[0]) * ry_to_ev

                    for key, value in [['energy', En], ['energy_accuracy', E_acc]]:
                        trajectory_data[key].append(value)
                        parsed_data[key + units_suffix] = default_energy_units
                    # TODO: decide units for magnetization. now bohr mag/cell
                    j = 0
                    while True:
                        j += 1
                        line2 = data_step[count + j]

                        for string, key in [
                                ['one-electron contribution', 'energy_one_electron'],
                                ['hartree contribution', 'energy_hartree'],
                                ['xc contribution', 'energy_xc'],
                                ['ewald contribution', 'energy_ewald'],
                                ['smearing contrib.', 'energy_smearing'],
                                ['one-center paw contrib.', 'energy_one_center_paw'],
                                ['est. exchange err', 'energy_est_exchange'],
                                ['Fock energy', 'energy_fock'],
                                # Add also ENVIRON specific contribution to the total energy
                                ['solvation energy', 'energy_solvation'],
                                ['cavitation energy', 'energy_cavitation'],
                                ['PV energy', 'energy_pv'],
                                ['periodic energy correct.', 'energy_pbc_correction'],
                                ['ionic charge energy', 'energy_ionic_charge'],
                                ['external charges energy', 'energy_external_charges']]:
                            if string in line2:
                                value = grep_energy_from_line(line2)
                                try:
                                    trajectory_data[key].append(value)
                                except KeyError:
                                    trajectory_data[key] = [value]
                                parsed_data[key + units_suffix] = default_energy_units
                        # magnetizations
                        if 'total magnetization' in line2:
                            this_m = line2.split('=')[1].split('Bohr')[0]
                            try:  # magnetization might be a scalar
                                value = float(this_m)
                            except ValueError:  # but can also be a three vector component in non-collinear calcs
                                value = [float(i) for i in this_m.split()]
                            try:
                                trajectory_data['total_magnetization'].append(value)
                            except KeyError:
                                trajectory_data['total_magnetization'] = [value]
                            parsed_data['total_magnetization' + units_suffix] = default_magnetization_units
                        elif 'absolute magnetization' in line2:
                            value = float(line2.split('=')[1].split('Bohr')[0])
                            try:
                                trajectory_data['absolute_magnetization'].append(value)
                            except KeyError:
                                trajectory_data['absolute_magnetization'] = [value]
                            parsed_data['absolute_magnetization' + units_suffix] = default_magnetization_units
                        # exit loop
                        elif 'convergence' in line2:
                            break

                    if vdw_correction:
                        j = 0
                        while True:
                            j += -1
                            line2 = data_step[count + j]
                            if 'Non-local correlation energy' in line2:
                                value = grep_energy_from_line(line2)
                                try:
                                    trajectory_data['energy_vdw'].append(value)
                                except KeyError:
                                    trajectory_data['energy_vdw'] = [value]
                                break
                        parsed_data['energy_vdw' + units_suffix] = default_energy_units
                except Exception:
                    logs.warning.append('Error while parsing for energy terms.')

            elif 'the Fermi energy is' in line:
                try:
                    value = float(line.split('is')[1].split('ev')[0])
                    try:
                        trajectory_data['fermi_energy'].append(value)
                    except KeyError:
                        trajectory_data['fermi_energy'] = [value]
                    parsed_data['fermi_energy' + units_suffix] = default_energy_units
                except Exception:
                    logs.warning.append('Error while parsing Fermi energy from the output file.')

            elif 'Forces acting on atoms' in line:
                try:
                    forces = []
                    j = 0
                    while True:
                        j += 1
                        line2 = data_step[count + j]
                        if 'atom ' in line2:
                            line2 = line2.split('=')[1].split()
                            # CONVERT FORCES IN eV/Ang
                            vec = [float(s) * ry_to_ev / bohr_to_ang for s in line2]
                            forces.append(vec)
                        if len(forces) == nat:
                            break
                    try:
                        trajectory_data['forces'].append(forces)
                    except KeyError:
                        trajectory_data['forces'] = [forces]
                    parsed_data['forces' + units_suffix] = default_force_units
                except Exception:
                    logs.warning.append('Error while parsing forces.')

            # TODO: adding the parsing support for the decomposition of the forces

            elif 'Total force =' in line:
                try:  # note that I can't check the units: not written in output!
                    value = float(line.split('=')[1].split('Total')[0]) * ry_to_ev / bohr_to_ang
                    try:
                        trajectory_data['total_force'].append(value)
                    except KeyError:
                        trajectory_data['total_force'] = [value]
                    parsed_data['total_force' + units_suffix] = default_force_units
                except Exception:
                    logs.warning.append('Error while parsing total force.')

            elif ('entering subroutine stress ...' in line) or ('Computing stress (Cartesian axis) and pressure' in line):
                try:
                    stress = []
                    for k in range(10 + 5 * vdw_correction):
                        if "P=" in data_step[count + k + 1]:
                            count2 = count + k + 1
                    if '(Ry/bohr**3)' not in data_step[count2]:
                        raise QEOutputParsingError('Error while parsing stress: unexpected units.')
                    for k in range(3):
                        line2 = data_step[count2 + k + 1].split()
                        vec = [float(s) * 10**(-9) * ry_si / (bohr_si)**3 for s in line2[0:3]]
                        stress.append(vec)
                    try:
                        trajectory_data['stress'].append(stress)
                    except KeyError:
                        trajectory_data['stress'] = [stress]
                    parsed_data['stress' + units_suffix] = default_stress_units
                except Exception:
                    logs.warning.append('Error while parsing stress tensor.')

            # Electronic and ionic dipoles when 'lelfield' was set to True in input parameters
            elif lelfield is True:

                if 'Electronic Dipole per cell' in line:
                    electronic_dipole = float(line.split()[-1])
                    trajectory_frame.setdefault('electronic_dipole_cell_average', []).append(electronic_dipole)

                elif 'Ionic Dipole per cell' in line:
                    ionic_dipole = float(line.split()[-1])
                    trajectory_frame.setdefault('ionic_dipole_cell_average', []).append(ionic_dipole)

                elif 'Electronic Dipole on Cartesian axes' in line:
                    electronic_dipole = [float(data_step[count + i + 1].split()[1]) for i in range(3)]
                    trajectory_frame.setdefault('electronic_dipole_cartesian_axes', []).append(electronic_dipole)

                elif 'Ionic Dipole on Cartesian axes' in line:
                    ionic_dipole = [float(data_step[count + i + 1].split()[1]) for i in range(3)]
                    trajectory_frame.setdefault('ionic_dipole_cartesian_axes', []).append(ionic_dipole)

        # End of trajectory frame, only keep last entries for dipole related values
        if lelfield is True:

            # For every property only get the last entry if possible
            try:
                ed_cell = trajectory_frame['electronic_dipole_cell_average'].pop()
            except IndexError:
                ed_cell = None

            try:
                ed_axes = trajectory_frame['electronic_dipole_cartesian_axes'].pop()
            except IndexError:
                ed_axes = None

            try:
                id_cell = trajectory_frame['ionic_dipole_cell_average'].pop()
            except IndexError:
                id_cell = None

            try:
                id_axes = trajectory_frame['ionic_dipole_cartesian_axes'].pop()
            except IndexError:
                id_axes = None

            # Only add them if all four properties were successfully parsed
            if all([value is not None for value in [ed_cell, ed_axes, id_cell, id_axes]]):
                trajectory_data.setdefault('electronic_dipole_cell_average', []).append(ed_cell)
                trajectory_data.setdefault('electronic_dipole_cartesian_axes', []).append(ed_axes)
                trajectory_data.setdefault('ionic_dipole_cell_average', []).append(id_cell)
                trajectory_data.setdefault('ionic_dipole_cartesian_axes', []).append(id_axes)

    # If specified in the parser options, parse the atomic occupations
    parse_atomic_occupations = parser_options.get('parse_atomic_occupations', False)

    if parse_atomic_occupations:

        atomic_occupations = {}
        hubbard_blocks = stdout.split('LDA+U parameters')

        for line in hubbard_blocks[-1].split('\n'):

            if 'Tr[ns(na)]' in line:

                values = line.split('=')
                atomic_index = values[0].split()[1]
                occupations = values[1].split()

                if len(occupations) == 1:
                    atomic_occupations[atomic_index] = {
                        'total': occupations[0]
                    }
                elif len(occupations) == 3:
                    atomic_occupations[atomic_index] = {
                        'up': occupations[0],
                        'down': occupations[1],
                        'total': occupations[2]
                    }
                else:
                    continue

        parsed_data['atomic_occupations'] = atomic_occupations

    # Ionic calculations and BFGS algorithm did not print that calculation is converged
    if 'atomic_positions_relax' in trajectory_data and not marker_bfgs_converged:
        logs.error.append('ERROR_IONIC_CONVERGENCE_NOT_REACHED')

    # Ionic calculation that hit the maximum number of ionic steps. Note: does not necessarily mean that convergence was
    # not reached as it could have occurred in the last step.
    if maximum_ionic_steps is not None and maximum_ionic_steps == parsed_data.get('number_ionic_steps', None):
        logs.warning.append('ERROR_MAXIMUM_IONIC_STEPS_REACHED')

    # Remove duplicate log messages by turning it into a set. Then convert back to list as that is what is expected
    logs.error = list(set(logs.error))
    logs.warning = list(set(logs.warning))

    parsed_data['bands'] = bands_data
    parsed_data['structure'] = structure_data
    parsed_data['trajectory'] = trajectory_data

    return parsed_data, logs


def grep_energy_from_line(line):
    try:
        return float(line.split('=')[1].split('Ry')[0]) * ry_to_ev
    except Exception:
        raise QEOutputParsingError('Error while parsing energy')


def convert_qe_time_to_sec(timestr):
    """Given the walltime string of Quantum Espresso, converts it in a number of seconds (float)."""
    rest = timestr.strip()

    if 'd' in rest:
        days, rest = rest.split('d')
    else:
        days = '0'

    if 'h' in rest:
        hours, rest = rest.split('h')
    else:
        hours = '0'

    if 'm' in rest:
        minutes, rest = rest.split('m')
    else:
        minutes = '0'

    if 's' in rest:
        seconds, rest = rest.split('s')
    else:
        seconds = '0.'

    if rest.strip():
        raise ValueError("Something remained at the end of the string '{}': '{}'".format(timestr, rest))

    num_seconds = (float(seconds) + float(minutes) * 60. + float(hours) * 3600. + float(days) * 86400.)

    return num_seconds


def parse_QE_errors(lines, line_number_start, warnings):
    """
    Parse QE errors messages (those appearing between some lines with ``%%%%%%%%``)

    :param lines: list of strings, the output text file as read by ``readlines()``
        or as obtained by ``data.split('\\\\n')`` when data is the text file read by ``read()``
    :param line_number_start: the line at which we identified some ``%%%%%%%%``
    :param warnings: dictionary where keys are error markers and the value the corresponding warning messages that
        should be returned.
    :return messages: a list of QE error messages
    """
    messages = []

    for line_number, line in enumerate(lines[line_number_start + 1:]):
        if '%%%%%%%%%%%%' in line:
            line_number_end = line_number
            break
    else:
        return messages

    for line in lines[line_number_start:line_number_end]:
        for marker, message in warnings.items():
            if marker in line:
                messages.append(message)

    return set(messages)
