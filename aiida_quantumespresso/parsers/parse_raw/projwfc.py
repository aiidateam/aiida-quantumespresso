from __future__ import absolute_import
import re


def parse_lowdin_charges(out_file_lines):
    """Parse the Lowdin Charge section of the ``projwfc.x`` output.

    :param out_file_lines: The projwfc.x stdout file content
    :type out_file_lines: list[str]
    :return: A tuple of lowdin data dict and spill parameter float
    """
    # find start of Lowdin charge output
    start_line = None
    for i, line in enumerate(out_file_lines):
        if line.strip() == 'Lowdin Charges:':
            if start_line is not None:
                raise IOError("'Lowdin Charges:' found on multiple lines: {}, {}".format(start_line, i))
            start_line = i

    if start_line is None:
        return None, None

    spill_parameter = None
    started_section = False
    atom_index = None
    output = {}
    for i, line in enumerate(out_file_lines[start_line + 1:]):
        lineno = start_line + 2 + i
        # break on empty line or spill parameter
        if not line.strip():
            if started_section:
                break
            started_section = True
            continue
        if line.strip().startswith('Spilling Parameter'):
            spill_parameter = float(line.strip().split()[-1])
            break
        atom_line_match = re.match(r'^\s*Atom\s\#\s*(\d+)\:(.*)$', line)
        if atom_line_match:
            atom_index = int(atom_line_match.group(1))
            line = atom_line_match.group(2)
        if atom_index is None:
            raise IOError('No atom index specified on or before line {}'.format(lineno))
        charge_type_match = re.match(r'^\s*(total charge|spin up|spin down|polarization)\s*\=\s*([\d\.]+)', line)
        if charge_type_match:
            charge_type = charge_type_match.group(1)
            output.setdefault(atom_index, {}).setdefault('sum', {})[charge_type] = float(charge_type_match.group(2))
        else:
            raise IOError('No charge type specified on line {}'.format(lineno))
        for orbital_type, value in re.findall(
            r'(s|p|d|f|px|py|pz|dxy|dxz|dyz|dz2|dx2-y2|f[xyz23\-\+]+)\s*\=\s*([\d\.]+)', line
        ):
            output.setdefault(atom_index, {}).setdefault(charge_type, {})[orbital_type] = float(value)

    return output, spill_parameter
