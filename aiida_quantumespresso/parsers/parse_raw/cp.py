# -*- coding: utf-8 -*-
from xml.dom.minidom import parseString

from xmlschema.etree import ElementTree

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_xml.cp.legacy import parse_cp_xml_output
from aiida_quantumespresso.parsers.parse_xml.parse import parse_xml_post_6_2
from aiida_quantumespresso.parsers.parse_xml.pw.legacy import parse_xml_child_integer
from aiida_quantumespresso.parsers.parse_xml.versions import QeXmlVersion, get_xml_file_version


def parse_cp_traj_stanzas(num_elements, splitlines, prepend_name, rescale=1.):
    """
    num_elements: Number of lines (with three elements) between lines with two only
    elements (containing step number and time in ps).
    num_elements is 3 for cell, and the number of atoms for coordinates and positions.

    splitlines: a list of lines of the file, already split in pieces using string.split

    prepend_name: a string to be prepended to the name of keys returned
    in the return dictionary.

    rescale: the values in each stanza are multiplied by this factor, for units conversion
    """
    steps = []
    times = []
    stanzas = []
    this_stanza = []
    start_stanza = False
    linenum = -1
    try:
        for linenum, l in enumerate(splitlines):
            if len(l) == 2:
                steps.append(int(l[0]))
                if set(l[1]) == {'*'}:
                    times.append(-1.0)
                else:
                    times.append(float(l[1]))
                start_stanza = True
                if len(this_stanza) != 0:
                    raise ValueError('Wrong position of short line.')
            elif len(l) == 3:
                if len(this_stanza) == 0 and not start_stanza:
                    raise ValueError('Wrong position of long line.')
                start_stanza = False
                this_stanza.append([float(l[0]) * rescale, float(l[1]) * rescale, float(l[2]) * rescale])
                if len(this_stanza) == num_elements:
                    stanzas.append(this_stanza)
                    this_stanza = []
            else:
                raise ValueError(f'Wrong line length ({len(l)})')
        if len(this_stanza) != 0:
            raise ValueError(f'Wrong length of last block ({len(this_stanza)} lines instead of 0).')
        if len(steps) != len(stanzas):
            raise ValueError('Length mismatch between number of steps and number of defined stanzas.')
        return {
            f'{prepend_name}_steps': steps,
            f'{prepend_name}_times': times,
            f'{prepend_name}_data': stanzas,
        }
    except Exception as e:
        e.message = f'At line {linenum + 1}: {e}'
        raise e


def parse_cp_text_output(data, xml_data):
    """data must be a list of strings, one for each lines, as returned by readlines().

    On output, a dictionary with parsed values
    """
    # TODO: uniform readlines() and read() usage for passing input to the parser

    parsed_data = {}
    parsed_data['warnings'] = []
    conjugate_gradient = False
    for count, line in enumerate(data):
        if 'conjugate gradient' in line.lower():
            conjugate_gradient = True
            parsed_data['conjugate_gradient'] = True
        if 'warning' in line.lower():
            parsed_data['warnings'].append(line)
        elif 'bananas' in line:
            parsed_data['warnings'].append('Bananas from the ortho.')
        elif 'CP' in line and 'WALL' in line:
            try:
                time = line.split('CPU')[1].split('WALL')[0]
                parsed_data['wall_time'] = time
            except:
                raise QEOutputParsingError('Error while parsing wall time.')
    #when the cp does a cg, the output is different and the parser below does not work
    #TODO: understand what the cg prints out and parse it (it is undocumented)
    if not conjugate_gradient:
        for count, line in enumerate(reversed(data)):
            if 'nfi' in line and 'ekinc' in line and 'econs' in line:
                this_line = data[len(data) - count]
                try:
                    parsed_data['ekinc'] = [float(this_line.split()[1])]
                except ValueError:
                    pass
                try:
                    parsed_data['temph'] = [float(this_line.split()[2])]
                except ValueError:
                    pass
                try:
                    parsed_data['tempp'] = [float(this_line.split()[3])]
                except ValueError:
                    pass
                try:
                    parsed_data['etot'] = [float(this_line.split()[4])]
                except ValueError:
                    pass
                try:
                    parsed_data['enthal'] = [float(this_line.split()[5])]
                except ValueError:
                    pass
                try:
                    parsed_data['econs'] = [float(this_line.split()[6])]
                except ValueError:
                    pass
                try:
                    parsed_data['econt'] = [float(this_line.split()[7])]
                except ValueError:
                    pass
                try:
                    parsed_data['vnhh'] = [float(this_line.split()[8])]
                except (ValueError, IndexError):
                    pass
                try:
                    parsed_data['xnhh0'] = [float(this_line.split()[9])]
                except (ValueError, IndexError):
                    pass
                try:
                    parsed_data['vnhp'] = [float(this_line.split()[10])]
                except (ValueError, IndexError):
                    pass
                try:
                    parsed_data['xnhp0'] = [float(this_line.split()[11])]
                except (ValueError, IndexError):
                    pass

    return parsed_data


def parse_cp_xml_counter_output(data):
    """Parse xml file print_counter.xml data must be a single string, as returned by file.read() (notice the difference
    with parse_text_output!) On output, a dictionary with parsed values."""
    dom = parseString(data)
    parsed_data = {}
    cardname = 'LAST_SUCCESSFUL_PRINTOUT'

    card1 = [_ for _ in dom.childNodes if _.nodeName == 'PRINT_COUNTER'][0]
    card2 = [_ for _ in card1.childNodes if _.nodeName == 'LAST_SUCCESSFUL_PRINTOUT'][0]

    tagname = 'STEP'
    parsed_data[cardname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, card2)

    return parsed_data


def parse_cp_counter_output(data):
    """Parse file print_counter data must be a single string, as returned by file.read() (notice the difference
    with parse_text_output!) On output, a dictionary with parsed values."""
    parsed_data = {}
    cardname = 'LAST_SUCCESSFUL_PRINTOUT'
    tagname = 'STEP'
    numbers = [int(s) for s in data.split() if s.isdigit()]
    if numbers:
        parsed_data[cardname.lower().replace('-', '_')] = numbers[0]

    return parsed_data


def parse_cp_raw_output(stdout, output_xml=None, xml_counter_file=None, print_counter_xml=True):

    parser_warnings = []

    # analyze the xml
    if output_xml is not None:
        xml_parsed = ElementTree.ElementTree(element=ElementTree.fromstring(output_xml))
        xml_file_version = get_xml_file_version(xml_parsed)
        if xml_file_version == QeXmlVersion.POST_6_2:
            xml_data, logs = parse_xml_post_6_2(xml_parsed)
        elif xml_file_version == QeXmlVersion.PRE_6_2:
            xml_data = parse_cp_xml_output(output_xml)
    else:
        parser_warnings.append('Skipping the parsing of the xml file.')
        xml_data = {}

    # analyze the counter file, which keeps info on the steps
    if xml_counter_file is not None:
        if print_counter_xml:
            xml_counter_data = parse_cp_xml_counter_output(xml_counter_file)
        else:
            xml_counter_data = parse_cp_counter_output(xml_counter_file)
    else:
        xml_counter_data = {}
    stdout = stdout.split('\n')
    # understand if the job ended smoothly
    job_successful = any('JOB DONE' in line for line in reversed(stdout))

    out_data = parse_cp_text_output(stdout, xml_data)

    for key in out_data.keys():
        if key in list(xml_data.keys()):
            raise AssertionError(f'{key} found in both dictionaries')
        if key in list(xml_counter_data.keys()):
            raise AssertionError(f'{key} found in both dictionaries')
        # out_data keys take precedence and overwrite xml_data keys,
        # if the same key name is shared by both (but this should not happen!)

    final_data = {
        **xml_data,
        **out_data,
        **xml_counter_data,
    }

    if parser_warnings:
        final_data['parser_warnings'] = parser_warnings

    # TODO: parse the trajectory and save them in a reasonable format

    return final_data, job_successful
