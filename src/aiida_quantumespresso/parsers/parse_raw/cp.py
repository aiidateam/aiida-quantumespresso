import contextlib
from xml.dom.minidom import parseString
from xml.etree import ElementTree

from aiida_quantumespresso.parsers.parse_xml.cp.legacy import parse_cp_xml_output
from aiida_quantumespresso.parsers.parse_xml.parse import parse_xml_post_6_2
from aiida_quantumespresso.parsers.parse_xml.pw.legacy import parse_xml_child_integer
from aiida_quantumespresso.parsers.parse_xml.versions import QeXmlVersion, get_xml_file_version


def parse_cp_traj_stanzas(num_elements, splitlines, prepend_name, rescale=1.0):
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
        for index, line in enumerate(splitlines):
            linenum = index
            if len(line) == 2:
                steps.append(int(line[0]))
                if set(line[1]) == {'*'}:
                    times.append(-1.0)
                else:
                    times.append(float(line[1]))
                start_stanza = True
                if len(this_stanza) != 0:
                    raise ValueError('Wrong position of short line.')
            elif len(line) == 3:
                if len(this_stanza) == 0 and not start_stanza:
                    raise ValueError('Wrong position of long line.')
                start_stanza = False
                this_stanza.append([float(line[0]) * rescale, float(line[1]) * rescale, float(line[2]) * rescale])
                if len(this_stanza) == num_elements:
                    stanzas.append(this_stanza)
                    this_stanza = []
            else:
                raise ValueError(f'Wrong line length ({len(line)})')
        if len(this_stanza) != 0:
            raise ValueError(f'Wrong length of last block ({len(this_stanza)} lines instead of 0).')
        if len(steps) != len(stanzas):
            raise ValueError('Length mismatch between number of steps and number of defined stanzas.')
    except Exception as e:
        e.message = f'At line {linenum + 1}: {e}'
        raise

    return {
        f'{prepend_name}_steps': steps,
        f'{prepend_name}_times': times,
        f'{prepend_name}_data': stanzas,
    }


def parse_cp_text_output(stdout):
    """stdout must be a list of strings, one for each lines, as returned by readlines().

    On output, a dictionary with parsed values
    """
    # TODO: uniform readlines() and read() usage for passing input to the parser

    parsed_data = {}
    parsed_data['warnings'] = []
    conjugate_gradient = False
    for line in stdout:
        if 'conjugate gradient' in line.lower():
            conjugate_gradient = True
            parsed_data['conjugate_gradient'] = True
        if 'warning' in line.lower():
            parsed_data['warnings'].append(line)
        elif 'bananas' in line:
            parsed_data['warnings'].append('Bananas from the ortho.')

    # when the cp does a cg, the output is different and the parser below does not work
    # TODO: understand what the cg prints out and parse it (it is undocumented)
    if not conjugate_gradient:
        for count, line in enumerate(reversed(stdout)):
            if 'nfi' in line and 'ekinc' in line and 'econs' in line:
                this_line = stdout[len(stdout) - count]
                with contextlib.suppress(ValueError):
                    parsed_data['ekinc'] = [float(this_line.split()[1])]
                with contextlib.suppress(ValueError):
                    parsed_data['temph'] = [float(this_line.split()[2])]
                with contextlib.suppress(ValueError):
                    parsed_data['tempp'] = [float(this_line.split()[3])]
                with contextlib.suppress(ValueError):
                    parsed_data['etot'] = [float(this_line.split()[4])]
                with contextlib.suppress(ValueError):
                    parsed_data['enthal'] = [float(this_line.split()[5])]
                with contextlib.suppress(ValueError):
                    parsed_data['econs'] = [float(this_line.split()[6])]
                with contextlib.suppress(ValueError):
                    parsed_data['econt'] = [float(this_line.split()[7])]
                with contextlib.suppress(ValueError, IndexError):
                    parsed_data['vnhh'] = [float(this_line.split()[8])]
                with contextlib.suppress(ValueError, IndexError):
                    parsed_data['xnhh0'] = [float(this_line.split()[9])]
                with contextlib.suppress(ValueError, IndexError):
                    parsed_data['vnhp'] = [float(this_line.split()[10])]
                with contextlib.suppress(ValueError, IndexError):
                    parsed_data['xnhp0'] = [float(this_line.split()[11])]

    return parsed_data


def parse_cp_xml_counter_output(data):
    """Parse xml file print_counter.xml data must be a single string, as returned by file.read() (notice the difference
    with parse_text_output!) On output, a dictionary with parsed values."""
    dom = parseString(data)
    parsed_data = {}
    cardname = 'LAST_SUCCESSFUL_PRINTOUT'

    card1 = next(_ for _ in dom.childNodes if _.nodeName == 'PRINT_COUNTER')
    card2 = next(_ for _ in card1.childNodes if _.nodeName == 'LAST_SUCCESSFUL_PRINTOUT')

    tagname = 'STEP'
    parsed_data[cardname.lower().replace('-', '_')] = parse_xml_child_integer(tagname, card2)

    return parsed_data


def parse_cp_counter_output(data):
    """Parse file print_counter data must be a single string, as returned by file.read() (notice the difference
    with parse_text_output!) On output, a dictionary with parsed values."""
    parsed_data = {}
    cardname = 'LAST_SUCCESSFUL_PRINTOUT'
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

    out_data = parse_cp_text_output(stdout)

    for key in out_data:
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
