# -*- coding: utf-8 -*-
"""A collection of function that are used to parse the output of Quantum Espresso PHonon.

The function that needs to be called from outside is parse_raw_ph_output(). Ideally, the functions should work even
without aiida and will return a dictionary with parsed keys.
"""
from __future__ import annotations

import numpy
from qe_tools import CONSTANTS

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.base import convert_qe_time_to_sec
from aiida_quantumespresso.parsers.parse_xml.pw.legacy import parse_xml_child_bool, read_xml_card


def parse_raw_ph_output(stdout, logs, tensors=None, dynamical_matrices=None):
    """Parses the raw output of a Quantum ESPRESSO `ph.x` calculation.

    :param stdout: the content of the stdout file as a string
    :param tensors: the content of the tensors.xml file as a string
    :param dynamical_matrices: a list of the content of the dynamical matrix files as a string
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    data_lines = stdout.split('\n')

    # Parse tensors, if present
    tensor_data = {}
    if tensors:
        try:
            tensor_data = parse_ph_tensor(tensors)
        except QEOutputParsingError:
            logs.warning.append('Error while parsing the tensor files')

    out_data = parse_ph_text_output(data_lines, logs)

    # parse dynamical matrices if present
    dynmat_data = {}
    if dynamical_matrices:
        # find lattice parameter
        for dynmat_counter, dynmat in enumerate(dynamical_matrices):

            lines = dynmat.split('\n')

            # check if the file contains frequencies (i.e. is useful) or not
            dynmat_to_parse = False
            if not lines:
                continue
            try:
                _ = [float(i) for i in lines[0].split()]
            except ValueError:
                dynmat_to_parse = True

            if not dynmat_to_parse:
                continue

            # parse it
            this_dynmat_data = parse_ph_dynmat(lines, logs)

            # join it with the previous dynmat info
            dynmat_data[f'dynamical_matrix_{dynmat_counter}'] = this_dynmat_data
            # TODO: use the bands format?

    # join dictionaries, there should not be any twice repeated key
    for key in out_data.keys():
        if key in list(tensor_data.keys()):
            raise AssertionError(f'{key} found in two dictionaries')
        if key in list(dynmat_data.keys()):
            raise AssertionError(f'{key} found in two dictionaries')

    symmetry_labels = out_data.pop('symmetry_labels', {})

    if len(dynmat_data) == 1 and 'dynamical_matrix_0' in dynmat_data:
        dynmat_data['dynamical_matrix_1'] = dynmat_data.pop('dynamical_matrix_0')

    parsed_data = dict(list(dynmat_data.items()) + list(out_data.items()) + list(tensor_data.items()))

    for q_index, q_symlabels in symmetry_labels.items():
        if f'dynamical_matrix_{q_index}' in parsed_data:
            parsed_data[f'dynamical_matrix_{q_index}'].update(q_symlabels)

    return parsed_data, logs


def parse_ph_tensor(data):
    """Parse the xml tensor file of QE v5.0.3 data must be read from the file with the .read() function (avoid
    readlines)"""
    from xml.dom.minidom import parseString

    dom = parseString(data)

    parsed_data = {}

    # card EF_TENSORS
    cardname = 'EF_TENSORS'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'DONE_ELECTRIC_FIELD'
    parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    if parsed_data[tagname.lower()]:
        try:
            second_tagname = 'DIELECTRIC_CONSTANT'
            parsed_data[second_tagname.lower()] = parse_xml_matrices(second_tagname, target_tags)
        except:
            raise QEOutputParsingError('Failed to parse Dielectric constant')

    tagname = 'DONE_EFFECTIVE_CHARGE_EU'
    parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    if parsed_data[tagname.lower()]:
        try:
            second_tagname = 'EFFECTIVE_CHARGES_EU'
            dumb_matrix = parse_xml_matrices(second_tagname, target_tags)
            # separate the elements of the messy matrix, with a matrix 3x3 for each element
            new_matrix = []
            this_at = []
            for i in dumb_matrix:
                this_at.append(i)
                if len(this_at) == 3:
                    new_matrix.append(this_at)
                    this_at = []

            parsed_data[second_tagname.lower()] = new_matrix
        except:
            raise QEOutputParsingError('Failed to parse effective charges eu')

    return parsed_data


def parse_xml_matrices(tagname, target_tags):
    """Can be used to load the disordered matrices of the QE XML file."""
    a = target_tags.getElementsByTagName(tagname)[0]
    b = a.childNodes[0]
    flat_array = b.data.split()
    # convert to float, then into a list of tuples, then into a list of lists
    flat_array = [float(i) for i in flat_array]
    list_tuple = list(zip(*[iter(flat_array)] * 3))
    return [list(i) for i in list_tuple]


def parse_ph_text_output(lines, logs):
    """Parses the stdout of Quantum ESPRESSO ``ph.x``.

    :param lines: list of strings, the file as read by readlines()
    :return: dictionary with parsed values
    """

    def parse_qpoints(lines):
        """Parse the q-points from the corresponding lines in the stdout."""

        return {int(line.split()[0]): [float(coord) for coord in line.split()[1:4]] for line in lines}

    def parse_mode_symmetries(lines, num_atoms):
        """Parse the mode symmetries from the block after diagonalization of the dynamical matrix."""

        q_data = {}
        symlabel_q_point = [float(i) for i in lines[2].split('q = (')[-1].split(')')[0].split()]

        q_data['mode_symmetry'] = []

        for line in lines:

            if 'Mode symmetry' in line:
                q_data['point_group'] = line.split('Mode symmetry,')[1].split('point group:')[0].strip()

            if '-->' in line:
                freq_start = int(line.split('(')[1].split(')')[0].split('-')[0])
                freq_end = int(line.split('(')[1].split(')')[0].split('-')[1])

                if line.split()[-1] == 'I':
                    symm_label = line.split()[-2]
                else:
                    symm_label = line.split()[-1]

                q_data['mode_symmetry'].extend([symm_label] * (freq_end - freq_start + 1))

        # If the mode symmetries are not printed, set those for the non-symmetry case
        if 'point_group' not in q_data:
            q_data['point_group'] = 'C_1'
            q_data['mode_symmetry'] = ['A'] * (3 * num_atoms)

        return symlabel_q_point, q_data

    parsed_data = {}

    parsed_data['num_q_found'] = 0

    # Parse number of q-points and number of atoms
    for count, line in enumerate(lines):

        if 'q-points for this run' in line:
            try:
                parsed_data['number_of_qpoints'] = int(line.split('/')[1].split('q-points')[0])
                parsed_data['q_points'] = parse_qpoints(lines[count + 2:count + parsed_data['number_of_qpoints'] + 2])

            except Exception as exc:
                logs.warning.append(f'Error while parsing number of q points: {exc}')

        elif 'q-points)' in line:
            # case of a 'only_wfc' calculation
            try:
                parsed_data['number_of_qpoints'] = int(line.split('q-points')[0].split('(')[1])
                parsed_data['q_points'] = parse_qpoints(lines[count + 2:count + parsed_data['number_of_qpoints'] + 2])

            except Exception as exc:
                logs.warning.append(f'Error while parsing number of q points: {exc}')

        # In case a single q-point is provided at the end of the input file, there is no summary of the q-points at the
        # start of the stdout. Then we have to parse the q-point further down the output file.
        elif 'number_of_qpoints' not in parsed_data and 'Calculation of q =' in line:
            parsed_data['number_of_qpoints'] = 1
            parsed_data['q_points'] = {1: [float(coord) for coord in line.split('=')[-1].split()]}

        elif 'number of atoms/cell' in line:
            try:
                num_atoms = int(line.split('=')[1])
                parsed_data['number_of_atoms'] = num_atoms
            except Exception:
                logs.warning.append('Error while parsing number of atoms.')

        elif 'irreducible representations' in line:
            if 'number_of_irr_representations_for_each_q' not in list(parsed_data.keys()):
                parsed_data['number_of_irr_representations_for_each_q'] = []
            try:
                num_irr_repr = int(line.split('irreducible')[0].split('are')[1])
                parsed_data['number_of_irr_representations_for_each_q'].append(num_irr_repr)
            except Exception:
                pass

        elif 'Diagonalizing the dynamical matrix' in line:

            mode_count = count

            while not 'Calculation' in lines[mode_count] and not 'CPU' in lines[mode_count]:
                mode_count += 1

            symlabel_q_point, q_data = parse_mode_symmetries(lines[count: mode_count], num_atoms)

            for q_index, q_point in parsed_data['q_points'].items():
                if numpy.isclose(
                    numpy.array(q_point), numpy.array(symlabel_q_point), rtol=0, atol=1e-7
                ).all():
                    parsed_data.setdefault('symmetry_labels', {})
                    parsed_data['symmetry_labels'][q_index] = q_data

            parsed_data['num_q_found'] += 1

    # Remove the q-points from the parsed data; these are only used to assign the symmetry labels to the right index
    parsed_data.pop('q_points', None)

    # Trim the number of irreps to the number of q-points finished in this run
    parsed_data['number_of_irr_representations_for_each_q'] = parsed_data.get(
        'number_of_irr_representations_for_each_q', []
    )[:parsed_data.pop('num_q_found')]

    return parsed_data


def parse_ph_dynmat(data, logs, lattice_parameter=None, also_eigenvectors=False, parse_header=False):
    """Parse frequencies and eigenvectors of a single dynamical matrix.

    :param data: the text read with the function readlines()
    :param lattice_parameter: the lattice_parameter ('alat' in QE jargon). If
        None, q_point is kept in 2pi/a coordinates as in the dynmat file.
    :param also_eigenvectors: if True, return an additional 'eigenvectors'
        array in output, containing also the eigenvectors. This will be
        a list of lists, that when converted to a numpy array has 4 indices,
        with shape Neigenstates x Natoms x 3(xyz) x 2 (re,im)
        To convert to a complex numpy array, you can use::

          ev = np.array(parsed_data['eigenvectors'])
          ev = ev[:,:,:,0] + 1j * ev[:,:,:,1]
    :param parse_header: if True, return additional keys in the returned
        parsed_data dictionary, including information from the header

    :return: a dictionary with parsed values and units
    """
    parsed_data = {}

    if 'Dynamical matrix file' not in data[0]:
        raise QEOutputParsingError('Dynamical matrix is not in the expected format')

    frequencies = []
    eigenvectors = []

    starting_line = 1
    if parse_header:
        header_dict = {'warnings': []}
        try:
            pieces = data[2].split()
            if len(pieces) != 9:
                raise QEOutputParsingError('Wrong # of elements on line 3')
            try:
                num_species = int(pieces[0])
                num_atoms = int(pieces[1])
                header_dict['ibrav'] = int(pieces[2])
                header_dict['celldm'] = [float(i) for i in pieces[3:]]
                # In angstrom
                alat = header_dict['celldm'][0] * CONSTANTS.bohr_to_ang
                if abs(alat) < 1.e-5:
                    raise QEOutputParsingError(
                        'Lattice constant=0! Probably you are using an '
                        'old Quantum ESPRESSO version?'
                    )
                header_dict['alat'] = alat
                header_dict['alat_units'] = 'angstrom'
            except ValueError:
                raise QEOutputParsingError('Wrong data on line 3')

            starting_line = 3
            if header_dict['ibrav'] == 0:
                if 'Basis vectors' not in data[3]:
                    raise QEOutputParsingError("Wrong format (no 'Basis vectors' line)")
                try:
                    v1 = [float(_) * alat for _ in data[4].split()]
                    v2 = [float(_) * alat for _ in data[5].split()]
                    v3 = [float(_) * alat for _ in data[6].split()]
                    if len(v1) != 3 or len(v2) != 3 or len(v3) != 3:
                        raise QEOutputParsingError('Wrong length for basis vectors')
                    header_dict['lattice_vectors'] = [v1, v2, v3]
                    header_dict['lattice_vectors_units'] = 'angstrom'
                except ValueError:
                    raise QEOutputParsingError('Wrong data for basis vectors')
                starting_line += 4

            species_info = {}
            species = []
            for idx, sp_line in enumerate(data[starting_line:starting_line + num_species], start=1):
                pieces = sp_line.split("'")
                if len(pieces) != 3:
                    raise QEOutputParsingError('Wrong # of elements for one of the species')
                try:
                    if int(pieces[0]) != idx:
                        raise QEOutputParsingError('Error with the indices of the species')
                    species.append([pieces[1].strip(), float(pieces[2]) / CONSTANTS.amu_Ry])
                except ValueError:
                    raise QEOutputParsingError('Error parsing the species')

            masses = dict(species)
            header_dict['masses'] = masses

            atoms_coords = []
            atoms_labels = []
            starting_line += num_species
            for idx, atom_line in enumerate(data[starting_line:starting_line + num_atoms], start=1):
                pieces = atom_line.split()
                if len(pieces) != 5:
                    raise QEOutputParsingError(
                        f'Wrong # of elements for one of the atoms: {len(pieces)}, line {starting_line + idx}: {pieces}'
                    )
                try:
                    if int(pieces[0]) != idx:
                        raise QEOutputParsingError(f'Error with the indices of the atoms: {int(pieces[0])} vs {idx}')
                    sp_idx = int(pieces[1])
                    if sp_idx > len(species):
                        raise QEOutputParsingError(f'Wrong index for the species: {sp_idx}, but max={len(species)}')
                    atoms_labels.append(species[sp_idx - 1][0])
                    atoms_coords.append([float(pieces[2]) * alat, float(pieces[3]) * alat, float(pieces[4]) * alat])
                except ValueError:
                    raise QEOutputParsingError('Error parsing the atoms')
                except IndexError:
                    raise QEOutputParsingError('Error with the indices in the atoms section')
            header_dict['atoms_labels'] = atoms_labels
            header_dict['atoms_coords'] = atoms_coords
            header_dict['atoms_coords_units'] = 'angstrom'

            starting_line += num_atoms

            starting_line += 1  # Got to the next line to check
            if 'Dynamical' not in data[starting_line]:
                raise QEOutputParsingError("Wrong format (no 'Dynamical  Matrix' line)")

            ## Here I finish the header parsing

        except QEOutputParsingError as e:
            logs.warning.append(
                'Problem parsing the header of the matdyn file! (msg: {}). '
                'Storing only the information I managed to retrieve'.format(e.message)
            )
            header_dict['warnings'].append(
                'There was some parsing error and this dictionary is '
                'not complete, see the warnings of the top parsed_data dict'
            )

        # I store what I got
        parsed_data['header'] = header_dict

    for line_counter, line in enumerate(data[starting_line:], start=starting_line):
        if 'q = ' in line:
            # q point is written several times, because it can also be rotated.
            # I consider only the first point, which is the one computed
            if 'q_point' not in parsed_data:
                q_point = [float(i) for i in line.split('(')[1].split(')')[0].split()]
                if lattice_parameter:
                    parsed_data['q_point'] = [e * 2 * numpy.pi / lattice_parameter for e in q_point]
                    parsed_data['q_point_units'] = 'angstrom-1'
                else:
                    parsed_data['q_point'] = q_point
                    parsed_data['q_point_units'] = '2pi/lattice_parameter'

        if 'freq' in line or 'omega' in line:
            this_freq = line.split('[cm-1]')[0].split('=')[-1]

            # exception for bad fortran coding: *** could be written instead of the number
            if '*' in this_freq:
                frequencies.append(None)
                logs.warning.append('Wrong fortran formatting found while parsing frequencies')
            else:
                frequencies.append(float(this_freq))

            this_eigenvectors = []
            for new_line in data[line_counter + 1:]:
                if (
                    'freq' in new_line or 'omega' in new_line or
                    '************************************************' in new_line
                ):
                    break
                this_things = new_line.split('(')[1].split(')')[0].split()
                try:
                    this_flatlist = [float(i) for i in this_things]
                except ValueError:
                    logs.warning.append('Wrong fortran formatting found while parsing eigenvectors')
                    # then save the three (xyz) complex numbers as [None,None]
                    this_eigenvectors.append([[None, None]] * 3)
                    continue

                list_tuples = list(zip(*[iter(this_flatlist)] * 2))
                # I save every complex number as a list of two numbers
                this_eigenvectors.append([[i[0], i[1]] for i in list_tuples])

            eigenvectors.append(this_eigenvectors)

    parsed_data['frequencies'] = frequencies
    parsed_data['frequencies_units'] = 'cm-1'
    # TODO: the eigenvectors should be written in the database according to a parser_opts.
    # for now, we don't store them, otherwise we get too much stuff
    # We implement anyway the possibility to get it with an optional parameter
    if also_eigenvectors:
        parsed_data['eigenvectors'] = eigenvectors

    return parsed_data


def parse_initialization_qpoints(stdout: str) -> dict:
    """Return the number of q-points from an initialization run.

    Here, the initialization run refers to the one performed by specifying
    `start_irr` and `last_irr` to 0 in the inputs.

    :return: parsed dictionary

    :raise: `RuntimeError` if the number of q-points cannot be parsed or it
        differs from the number of q-points in the stdout list.
    """
    import re

    parameters = {}

    # Regular expression to match `N` in `(  N q-points)`
    pattern = r'\(\s*(\d+)\s*q-points\)'
    match = re.search(pattern, stdout)
    if match:
        parameters.update({'number_of_qpoints': int(match.group(1))})
    else:
        raise RuntimeError('the number of q-points cannot be parsed')

    # Regular expression pattern to match the q-points section
    pattern = r'\(\s*\d+\s*q-points\):\s*\n\s*N\s*xq\(1\)\s*xq\(2\)\s*xq\(3\)\s*\n((?:\s*\d+\s*[\d\.\-\s]+\n?)*)'
    match = re.search(pattern, stdout)

    if match:
        q_points_block = match.group(1)

        # Regular expression to match each line of coordinates
        coord_pattern = r'\s*\d+\s*([\d\.\-]+)\s*([\d\.\-]+)\s*([\d\.\-]+)'

        coords = re.findall(coord_pattern, q_points_block) # Find all coordinates in the block
        q_points = [list(map(float, coord)) for coord in coords]
    else:
        raise RuntimeError('the list of q-points cannot be parsed')

    if parameters['number_of_qpoints'] != len(q_points):
        raise RuntimeError('the number of q-points do not coincde with the number of listed q-points')

    return parameters
