from __future__ import absolute_import
import re


def parse_lowdin_charges(out_file_lines, lowdin_lines=None):
    """Parse the Lowdin Charge section of the ``projwfc.x`` output.

    :param out_file_lines: The projwfc.x stdout file content
    :type out_file_lines: list[str]
    :param lowdin_lines: indices of lines where 'Lowdin Charges:' has been found
    :return: A tuple of lowdin data dict and spill parameter float
    """
    # find start of Lowdin charge output
    if lowdin_lines is None:
        lowdin_lines = []
        for i, line in enumerate(out_file_lines):
            if line.strip() == 'Lowdin Charges:':
                lowdin_lines.append(i)

    if len(lowdin_lines) > 1:
        raise IOError("'Lowdin Charges:' found on multiple lines: {}".format(lowdin_lines))

    if not lowdin_lines:
        return None, None

    spill_parameter = None
    started_section = False
    atom_index = None
    output = {}
    for i, line in enumerate(out_file_lines[lowdin_lines[0] + 1:]):
        lineno = lowdin_lines[0] + 2 + i
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
        charge_type_match = re.match(r'^\s*(total charge|spin up|spin down|polarization)\s*\=\s*(\-?[\d\.]+)', line)
        if charge_type_match:
            charge_type = charge_type_match.group(1).replace(' ', '_')
            output.setdefault(atom_index, {}).setdefault('sum', {})[charge_type] = float(charge_type_match.group(2))
        else:
            raise IOError('No charge type specified on line {}'.format(lineno))
        for orbital_type, value in re.findall(
            r'(s|p|d|f|px|py|pz|dxy|dxz|dyz|dz2|dx2-y2|f[xyz23\-\+]+)\s*\=\s*(\-?[\d\.]+)', line
        ):
            output.setdefault(atom_index, {}).setdefault(charge_type, {})[orbital_type] = float(value)

    return output, spill_parameter
