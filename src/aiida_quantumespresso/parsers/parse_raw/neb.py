# -*- coding: utf-8 -*-
"""A collection of function that are used to parse the output of Quantum Espresso Neb.

The function that needs to be called from outside is parse_raw_output_neb(). The functions mostly work without aiida
specific functionalities. The parsing will try to convert whatever it can in some dictionary, which by operative
decision doesn't have much structure encoded, [the values are simple ]
"""
from qe_tools import CONSTANTS

def detect_important_message(logs, line):
    message_map = {
        'error': {'invalid number of pools, out of range': 'ERROR_NPOOLS_TOO_HIGH',
                  'invalid number of images, out of range': 'ERROR_NIMAGE_HIGHER_THAN_NPROC',
                  'n. of images must be divisor of nproc': 'ERROR_NIMAGE_NOT_DIVISOR_OF_NPROC',
                  'nimage is larger than the available number of images': 'ERROR_NIMAGE_HIGHER_THAN_IMAGES'},
    }

    for marker, message in message_map['error'].items():
        if hasattr(marker, 'search'):
            if marker.match(line):
                if message is None:
                    message = line
                logs['errors'].append(message)
        else:
            if marker in line:
                if message is None:
                    message = line
                logs['errors'].append(message)
    
def parse_raw_output_neb(stdout):
    """Parses the output of a neb calculation Receives in input the paths to the output file.

    :param stdout: the stdout content as a string

    :return parameter_data: a dictionary with parsed parameters
    :return iteration_data: a dictionary with arrays (for relax & md calcs.)
    """
    import copy

    parser_warnings = []

    # parse the text output of the neb calculation
    out_data, iteration_data = parse_neb_text_output(stdout)

    # I add in the out_data all the last elements of iteration_data values.
    # I leave the possibility to skip some large arrays (None for the time being).
    skip_keys = []
    tmp_iteration_data = copy.copy(iteration_data)
    for k, v in tmp_iteration_data.items():
        if k in skip_keys:
            continue
        out_data[k] = v[-1]

    parameter_data = dict(list(out_data.items()) + [('parser_warnings', parser_warnings)])

    return parameter_data, iteration_data


def parse_neb_text_output(data):
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

    parsed_data = {}
    parsed_data['warnings'] = []
    parsed_data['errors'] = []
    iteration_data = defaultdict(list)

    # set by default the calculation as not converged.
    parsed_data['converged'] = [False, 0]

    for count, line in enumerate(data.split('\n')):
        detect_important_message(parsed_data, line)

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

    num_images = parsed_data['num_of_images'] if 'num_of_images' in parsed_data else None

    iteration_lines = data.split('-- iteration')[1:]
    iteration_lines = [i.split('\n') for i in iteration_lines]

    for iteration in iteration_lines:
        for count, line in enumerate(iteration):
            if 'activation energy (->)' in line:
                try:
                    activ_energy = float(line.split('=')[1].split('eV')[0])
                except Exception:
                    activ_energy = None
                iteration_data['forward_activation_energy'].append(activ_energy)
            elif 'activation energy (<-)' in line:
                try:
                    activ_energy = float(line.split('=')[1].split('eV')[0])
                except Exception:
                    activ_energy = None
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
                try:
                    path_length = float(line.split('=')[1].split('bohr')[0])
                    iteration_data['path_length'].append(path_length * CONSTANTS.bohr_to_ang)
                except Exception:
                    iteration_data['path_length'].append(None)
            elif 'inter-image distance' in line:
                try:
                    image_dist = float(line.split('=')[1].split('bohr')[0])
                    iteration_data['image_dist'].append(image_dist * CONSTANTS.bohr_to_ang)
                except Exception:
                    iteration_data['image_dist'].append(None)

    return parsed_data, dict(iteration_data)
