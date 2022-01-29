# -*- coding: utf-8 -*-
"""A collection of function that are used to parse the output of Quantum Espresso Neb.

The function that needs to be called from outside is parse_raw_output_neb(). The functions mostly work without aiida
specific functionalities. The parsing will try to convert whatever it can in some dictionary, which by operative
decision doesn't have much structure encoded, [the values are simple ]
"""
from qe_tools import CONSTANTS

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw import convert_qe_time_to_sec


def parse_raw_output_neb(stdout, input_dict, parser_opts=None):
    """Parses the output of a neb calculation Receives in input the paths to the output file.

    :param stdout: the stdout content as a string
    :param input_dict: dictionary with the neb input parameters
    :param parser_opts: not used

    :return parameter_data: a dictionary with parsed parameters
    :return iteration_data: a dictionary with arrays (for relax & md calcs.)
    :return job_successful: a boolean that is False in case of failed calculations

    :raises QEOutputParsingError: for errors in the parsing,

    2 different keys to check in output: parser_warnings and warnings.
    On an upper level, these flags MUST be checked.
    The first is expected to be empty unless QE failures or unfinished jobs.
    """
    import copy

    job_successful = True
    parser_warnings = []

    if not stdout:  # there is an output file, but it's empty -> crash
        job_successful = False

    # check if the job has finished (that doesn't mean without errors)
    finished_run = False
    for line in stdout.split('\n')[::-1]:
        if 'JOB DONE' in line:
            finished_run = True
            break
    if not finished_run:  # error if the job has not finished
        warning = 'QE neb run did not reach the end of the execution.'
        parser_warnings.append(warning)
        job_successful = False

    # parse the text output of the neb calculation
    try:
        out_data, iteration_data, critical_messages = parse_neb_text_output(stdout, input_dict)
    except QEOutputParsingError as exc:
        if not finished_run:  # I try to parse it as much as possible
            parser_warnings.append('Error while parsing the output file')
            out_data = {'warnings': []}
            iteration_data = {}
            critical_messages = []
        else:  # if it was finished and I got an error, it's a mistake of the parser
            raise QEOutputParsingError(f'Error while parsing NEB text output: {exc}')

    # I add in the out_data all the last elements of iteration_data values.
    # I leave the possibility to skip some large arrays (None for the time being).
    skip_keys = []
    tmp_iteration_data = copy.copy(iteration_data)
    for k, v in tmp_iteration_data.items():
        if k in skip_keys:
            continue
        out_data[k] = v[-1]

    # if there is a severe error, the calculation is FAILED
    if any([x in out_data['warnings'] for x in critical_messages]):
        job_successful = False

    parameter_data = dict(list(out_data.items()) + [('parser_warnings', parser_warnings)])

    # return various data.
    # parameter data will be mapped in Dict
    # iteration_data in ArrayData
    return parameter_data, iteration_data, job_successful


def parse_neb_text_output(data, input_dict={}):
    """Parses the text output of QE Neb.

    :param data: a string, the file as read by read()
    :param input_dict: dictionary with the input parameters

    :return parsed_data: dictionary with key values, referring to quantities
                         at the last step.
    :return iteration_data: key,values referring to intermediate iterations.
                             Empty dictionary if no value is present.
    :return critical_messages: a list with critical messages. If any is found in
                               parsed_data['warnings'], the calculation is FAILED!
    """
    from collections import defaultdict

    from aiida_quantumespresso.parsers.parse_raw import parse_output_error
    from aiida_quantumespresso.utils.mapping import get_logging_container

    # TODO: find a more exhaustive list of the common errors of neb
    # critical warnings: if any is found, the calculation status is FAILED
    critical_warnings = {
        'scf convergence NOT achieved on image': 'SCF did not converge for a given image',
        'Maximum CPU time exceeded': 'Maximum CPU time exceeded',
        'reached the maximum number of steps': 'Maximum number of iterations reached in the image optimization',
    }

    minor_warnings = {
        'Warning:': None,
    }

    all_warnings = dict(list(critical_warnings.items()) + list(minor_warnings.items()))

    parsed_data = {}
    parsed_data['warnings'] = []
    iteration_data = defaultdict(list)

    # parse time, starting from the end
    # apparently, the time is written multiple times
    for line in reversed(data.split('\n')):
        if 'NEB' in line and 'WALL' in line:
            try:
                time = line.split('CPU')[1].split('WALL')[0].strip()
                parsed_data['wall_time'] = time
            except Exception:
                parsed_data['warnings'].append('Error while parsing wall time.')

            try:
                parsed_data['wall_time_seconds'] = \
                    convert_qe_time_to_sec(parsed_data['wall_time'])
            except ValueError:
                raise QEOutputParsingError('Unable to convert wall_time in seconds.')
            break

    # set by default the calculation as not converged.
    parsed_data['converged'] = [False, 0]
    logs = get_logging_container()
    lines = data.split('\n')

    for count, line in enumerate(lines):
        if 'initial path length' in line:
            initial_path_length = float(line.split('=')[1].split('bohr')[0])
            parsed_data['initial_path_length'] = initial_path_length * CONSTANTS.bohr_to_ang
        elif 'initial inter-image distance' in line:
            initial_image_dist = float(line.split('=')[1].split('bohr')[0])
            parsed_data['initial_image_dist'] = initial_image_dist * CONSTANTS.bohr_to_ang
        elif 'string_method' in line:
            parsed_data['string_method'] = line.split('=')[1].strip()
        elif 'restart_mode' in line:
            parsed_data['restart_mode'] = line.split('=')[1].strip()
        elif 'opt_scheme' in line:
            parsed_data['opt_scheme'] = line.split('=')[1].strip()
        elif 'num_of_images' in line:
            parsed_data['num_of_images'] = int(line.split('=')[1])
        elif 'nstep_path' in line:
            parsed_data['nstep_path'] = int(line.split('=')[1])
        elif 'CI_scheme' in line:
            parsed_data['ci_scheme'] = line.split('=')[1].strip()
        elif 'first_last_opt' in line:
            parsed_data['first_last_opt'] = True if line.split('=')[1] == 'T' else False
        elif 'use_freezing' in line:
            parsed_data['use_freezing'] = True if line.split('=')[1] == 'T' else False
        elif ' ds ' in line:
            parsed_data['ds_au'] = float(line.split('=')[1].split('a.u.')[0])
        elif '   k_max' in line:
            parsed_data['k_max'] = float(line.split('=')[1].split('a.u.')[0])
        elif '   k_min_au' in line:
            parsed_data['k_min_au'] = float(line.split('=')[1].split('a.u.')[0])
        elif 'suggested k_max' in line:
            parsed_data['suggested_k_max_au'] = float(line.split('=')[1].split('a.u.')[0])
        elif 'suggested k_min' in line:
            parsed_data['suggested_k_min_au'] = float(line.split('=')[1].split('a.u.')[0])
        elif 'path_thr' in line:
            parsed_data['path_thr'] = float(line.split('=')[1].split('eV')[0])
        elif 'list of climbing images' in line:
            parsed_data['climbing_images_manual'] = [int(_) for _ in line.split(':')[1].split(',')[:-1]]
        elif 'neb: convergence achieved in' in line:
            parsed_data['converged'] = [True, int(line.split('iteration')[0].split()[-1])]
        elif '%%%%%%%%%%%%%%' in line:
            parse_output_error(lines, count, logs)
        elif any(i in line for i in all_warnings):
            message = [all_warnings[i] for i in all_warnings.keys() if i in line][0]

            if message is not None:
                parsed_data['warnings'].append(message)

    parsed_data['warnings'].extend(logs.error)

    try:
        num_images = parsed_data['num_of_images']
    except KeyError:
        try:
            num_images = input_dict['PATH']['num_of_images']
        except KeyError:
            raise QEOutputParsingError(
                'No information on the number '
                'of images available (neither in input nor in output'
            )

    iteration_lines = data.split('-- iteration')[1:]
    iteration_lines = [i.split('\n') for i in iteration_lines]

    for iteration in iteration_lines:
        for count, line in enumerate(iteration):
            if 'activation energy (->)' in line:
                activ_energy = float(line.split('=')[1].split('eV')[0])
                iteration_data['forward_activation_energy'].append(activ_energy)
            elif 'activation energy (<-)' in line:
                activ_energy = float(line.split('=')[1].split('eV')[0])
                iteration_data['backward_activation_energy'].append(activ_energy)
            elif 'image        energy (eV)        error (eV/A)        frozen' in line:
                energies = []
                forces = []
                frozen = []
                try:
                    for i in range(num_images):
                        split_line = iteration[count + 2 + i].split()[1:]
                        energies.append(float(split_line[0]))
                        forces.append(float(split_line[1]))
                        frozen.append(True if split_line[2] == 'T' else False)
                    iteration_data['image_energies'].append(energies)
                    iteration_data['image_forces'].append(forces)
                    iteration_data['image_frozen'].append(frozen)
                except Exception:
                    parsed_data['warnings'].append('Error while parsing the image energies and forces.')
            elif 'climbing image' in line:
                iteration_data['climbing_image_auto'].append([int(_) for _ in line.split('=')[1].split(',')])
            elif 'path length' in line:
                path_length = float(line.split('=')[1].split('bohr')[0])
                iteration_data['path_length'].append(path_length * CONSTANTS.bohr_to_ang)
            elif 'inter-image distance' in line:
                image_dist = float(line.split('=')[1].split('bohr')[0])
                iteration_data['image_dist'].append(image_dist * CONSTANTS.bohr_to_ang)

    return parsed_data, dict(iteration_data), list(critical_warnings.values())
