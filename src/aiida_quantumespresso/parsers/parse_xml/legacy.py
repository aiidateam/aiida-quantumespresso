# -*- coding: utf-8 -*-
"""Code that was written to parse the legacy XML format of Quantum ESPRESSO, which was deprecated in version 6.4."""
import string
from xml.dom.minidom import Element

from qe_tools import CONSTANTS

from aiida_quantumespresso.parsers import QEOutputParsingError

from .parse import cell_volume

units_suffix = '_units'
default_energy_units = 'eV'
default_k_points_units = '1 / angstrom'
default_length_units = 'Angstrom'


# In the following, some functions that helps the parsing of
# the xml file of QE v5.0.x (version below not tested)
def read_xml_card(dom, cardname):
    try:
        root_node = [_ for _ in dom.childNodes if isinstance(_, Element) and _.nodeName == 'Root'][0]
        the_card = [_ for _ in root_node.childNodes if _.nodeName == cardname][0]
        #the_card = dom.getElementsByTagName(cardname)[0]
        return the_card
    except Exception as e:
        print(e)
        raise QEOutputParsingError(f'Error parsing tag {cardname}')


def parse_xml_child_integer(tagname, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.childNodes[0]
        return int(b.data)
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}')


def parse_xml_child_float(tagname, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.childNodes[0]
        return float(b.data)
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}')


def parse_xml_child_bool(tagname, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.childNodes[0]
        return str2bool(b.data)
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}')


def str2bool(string):
    try:
        false_items = ['f', '0', 'false', 'no']
        true_items = ['t', '1', 'true', 'yes']
        string = str(string.lower().strip())
        if string in false_items:
            return False
        if string in true_items:
            return True
        else:
            raise QEOutputParsingError(f'Error converting string {string} to boolean value.')
    except Exception:
        raise QEOutputParsingError('Error converting string to boolean.')


def parse_xml_child_str(tagname, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.childNodes[0]
        return str(b.data).rstrip().replace('\n', '')
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}')


def parse_xml_child_attribute_str(tagname, attributename, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        value = str(a.getAttribute(attributename))
        return value.rstrip().replace('\n', '').lower()
    except Exception:
        raise QEOutputParsingError(
            f'Error parsing attribute {attributename}, tag {tagname} inside {target_tags.tagName}'
        )


def parse_xml_child_attribute_int(tagname, attributename, target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        value = int(a.getAttribute(attributename))
        return value
    except Exception:
        raise QEOutputParsingError(
            f'Error parsing attribute {attributename}, tag {tagname} inside {target_tags.tagName}'
        )


def convert_list_to_matrix(in_matrix, n_rows, n_columns):
    """converts a list into a list of lists (a matrix like) with n_rows and n_columns."""
    return [in_matrix[j:j + n_rows] for j in range(0, n_rows * n_columns, n_rows)]


def xml_card_cell(parsed_data, dom):
    #CARD CELL of QE output

    cardname = 'CELL'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['NON-PERIODIC_CELL_CORRECTION', 'BRAVAIS_LATTICE']:
        parsed_data[tagname.replace('-', '_').lower()] = parse_xml_child_str(tagname, target_tags)

    tagname = 'LATTICE_PARAMETER'
    value = parse_xml_child_float(tagname, target_tags)
    parsed_data[tagname.replace('-', '_').lower() + '_xml'] = value
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if metric not in ['bohr', 'angstrom']:
        raise QEOutputParsingError(
            f'Error parsing attribute {attrname}, tag {tagname} inside {target_tags.tagName}, units not found'
        )
    if metric == 'bohr':
        value *= CONSTANTS.bohr_to_ang
    parsed_data[tagname.replace('-', '_').lower()] = value

    tagname = 'CELL_DIMENSIONS'
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.childNodes[0]
        c = b.data.replace('\n', '').split()
        value = [float(i) for i in c]
        parsed_data[tagname.replace('-', '_').lower()] = value
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}.')

    tagname = 'DIRECT_LATTICE_VECTORS'
    lattice_vectors = []
    try:
        second_tagname = 'UNITS_FOR_DIRECT_LATTICE_VECTORS'
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b = a.getElementsByTagName('UNITS_FOR_DIRECT_LATTICE_VECTORS')[0]
        value = str(b.getAttribute('UNITS')).lower()
        parsed_data[second_tagname.replace('-', '_').lower()] = value

        metric = value
        if metric not in ['bohr', 'angstroms']:  # REMEMBER TO CHECK THE UNITS AT THE END OF THE FUNCTION
            raise QEOutputParsingError(
                f'Error parsing tag {tagname} inside {target_tags.tagName}: units not supported: {metric}'
            )

        lattice_vectors = []
        for second_tagname in ['a1', 'a2', 'a3']:
            #b = a.getElementsByTagName(second_tagname)[0]
            b = [_ for _ in a.childNodes if _.nodeName == second_tagname][0]
            c = b.childNodes[0]
            d = c.data.replace('\n', '').split()
            value = [float(i) for i in d]
            if metric == 'bohr':
                value = [CONSTANTS.bohr_to_ang * float(s) for s in value]
            lattice_vectors.append(value)

        volume = cell_volume(lattice_vectors[0], lattice_vectors[1], lattice_vectors[2])

    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName} inside {cardname}.')
    # NOTE: lattice_vectors will be saved later together with card IONS.atom

    tagname = 'RECIPROCAL_LATTICE_VECTORS'
    try:
        #a = target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

        second_tagname = 'UNITS_FOR_RECIPROCAL_LATTICE_VECTORS'
        b = a.getElementsByTagName(second_tagname)[0]
        value = str(b.getAttribute('UNITS')).lower()
        parsed_data[second_tagname.replace('-', '_').lower()] = value

        metric = value
        # NOTE: output is given in 2 pi / a [ang ^ -1]
        if metric not in ['2 pi / a']:
            raise QEOutputParsingError(
                f'Error parsing tag {tagname} inside {target_tags.tagName}: units {metric} not supported'
            )

        # reciprocal_lattice_vectors
        this_matrix = []
        for second_tagname in ['b1', 'b2', 'b3']:
            b = a.getElementsByTagName(second_tagname)[0]
            c = b.childNodes[0]
            d = c.data.replace('\n', '').split()
            value = [float(i) for i in d]
            if metric == '2 pi / a':
                value = [float(s) / parsed_data['lattice_parameter'] for s in value]
            this_matrix.append(value)
        parsed_data['reciprocal_lattice_vectors'] = this_matrix

    except Exception:
        raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}.')
    return parsed_data, lattice_vectors, volume


def xml_card_ions(parsed_data, dom, lattice_vectors, volume):
    cardname = 'IONS'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['NUMBER_OF_ATOMS', 'NUMBER_OF_SPECIES']:
        parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)

    tagname = 'UNITS_FOR_ATOMIC_MASSES'
    attrname = 'UNITS'
    parsed_data[tagname.lower()] = parse_xml_child_attribute_str(tagname, attrname, target_tags)

    try:
        parsed_data['species'] = {}
        parsed_data['species']['index'] = []
        parsed_data['species']['type'] = []
        parsed_data['species']['mass'] = []
        parsed_data['species']['pseudo'] = []
        for i in range(parsed_data['number_of_species']):
            tagname = 'SPECIE.' + str(i + 1)
            parsed_data['species']['index'].append(i + 1)

            #a=target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

            tagname2 = 'ATOM_TYPE'
            parsed_data['species']['type'].append(parse_xml_child_str(tagname2, a))

            tagname2 = 'MASS'
            parsed_data['species']['mass'].append(parse_xml_child_float(tagname2, a))

            tagname2 = 'PSEUDO'
            parsed_data['species']['pseudo'].append(parse_xml_child_str(tagname2, a))

        tagname = 'UNITS_FOR_ATOMIC_POSITIONS'
        attrname = 'UNITS'
        parsed_data[tagname.lower()] = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    except:
        raise QEOutputParsingError(f'Error parsing tag SPECIE.# inside {target_tags.tagName}.')


# TODO convert the units
# if parsed_data['units_for_atomic_positions'] not in ['alat','bohr','angstrom']:

    try:
        atomlist = []
        atoms_index_list = []
        atoms_if_pos_list = []
        tagslist = []
        for i in range(parsed_data['number_of_atoms']):
            tagname = 'ATOM.' + str(i + 1)
            # USELESS AT THE MOMENT, I DON'T SAVE IT
            # parsed_data['atoms']['list_index']=i
            #a=target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            tagname2 = 'INDEX'
            b = int(a.getAttribute(tagname2))
            atoms_index_list.append(b)
            tagname2 = 'SPECIES'

            chem_symbol = str(a.getAttribute(tagname2)).rstrip().replace('\n', '')
            # I check if it is a subspecie
            chem_symbol_digits = ''.join([i for i in chem_symbol if i in string.digits])
            try:
                tagslist.append(int(chem_symbol_digits))
            except ValueError:
                # If I can't parse the digit, it is probably not there: I add a None to the tagslist
                tagslist.append(None)
            # I remove the symbols
            chem_symbol = ''.join(i for i in chem_symbol if not i.isdigit())

            tagname2 = 'tau'
            b = a.getAttribute(tagname2)
            tau = [float(s) for s in b.rstrip().replace('\n', '').split()]
            metric = parsed_data['units_for_atomic_positions']
            if metric not in ['alat', 'bohr', 'angstrom']:  # REMEMBER TO CONVERT AT THE END
                raise QEOutputParsingError(f'Error parsing tag {tagname} inside {target_tags.tagName}')
            if metric == 'alat':
                tau = [parsed_data['lattice_parameter_xml'] * float(s) for s in tau]
            elif metric == 'bohr':
                tau = [CONSTANTS.bohr_to_ang * float(s) for s in tau]
            atomlist.append([chem_symbol, tau])
            tagname2 = 'if_pos'
            b = a.getAttribute(tagname2)
            if_pos = [int(s) for s in b.rstrip().replace('\n', '').split()]
            atoms_if_pos_list.append(if_pos)
        parsed_data['atoms'] = atomlist
        parsed_data['atoms_index_list'] = atoms_index_list
        parsed_data['atoms_if_pos_list'] = atoms_if_pos_list
        cell = {}
        cell['lattice_vectors'] = lattice_vectors
        cell['volume'] = volume
        cell['atoms'] = atomlist
        cell['tagslist'] = tagslist
        parsed_data['cell'] = cell
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag ATOM.# inside {target_tags.tagName}.')
    # saving data together with cell parameters. Did so for better compatibility with ASE.

    # correct some units that have been converted in
    parsed_data['atomic_positions' + units_suffix] = default_length_units
    parsed_data['direct_lattice_vectors' + units_suffix] = default_length_units

    return parsed_data


def xml_card_spin(parsed_data, dom):
    cardname = 'SPIN'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['LSDA', 'NON-COLINEAR_CALCULATION', 'SPIN-ORBIT_CALCULATION', 'SPIN-ORBIT_DOMAG']:
        parsed_data[tagname.replace('-', '_').lower()] = parse_xml_child_bool(tagname, target_tags)

    return parsed_data


def xml_card_header(parsed_data, dom):
    cardname = 'HEADER'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['FORMAT', 'CREATOR']:
        for attrname in ['NAME', 'VERSION']:
            parsed_data[(tagname + '_' + attrname).lower()
                        ] = parse_xml_child_attribute_str(tagname, attrname, target_tags)

    return parsed_data


def xml_card_planewaves(parsed_data, dom, calctype):
    if calctype not in ['pw', 'cp']:
        raise ValueError("Input flag not accepted, must be 'cp' or 'pw'")

    cardname = 'PLANE_WAVES'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'UNITS_FOR_CUTOFF'
    attrname = 'UNITS'
    units = parse_xml_child_attribute_str(tagname, attrname, target_tags).lower()
    if 'hartree' not in units:
        if 'rydberg' not in units:
            raise QEOutputParsingError(f'Units {units} are not supported by parser')
    else:
        if 'hartree' in units:
            conv_fac = CONSTANTS.hartree_to_ev
        else:
            conv_fac = CONSTANTS.ry_to_ev

        tagname = 'WFC_CUTOFF'
        parsed_data[tagname.lower()] = parse_xml_child_float(tagname, target_tags) * conv_fac
        parsed_data[tagname.lower() + units_suffix] = default_energy_units

        tagname = 'RHO_CUTOFF'
        parsed_data[tagname.lower()] = parse_xml_child_float(tagname, target_tags) * conv_fac
        parsed_data[tagname.lower() + units_suffix] = default_energy_units

    for tagname in ['FFT_GRID', 'SMOOTH_FFT_GRID']:
        grid = []
        for attrname in ['nr1', 'nr2', 'nr3']:
            if 'SMOOTH' in tagname:
                attrname += 's'
            grid.append(parse_xml_child_attribute_int(tagname, attrname, target_tags))
        parsed_data[tagname.lower()] = grid

    if calctype == 'cp':

        for tagname in ['MAX_NUMBER_OF_GK-VECTORS', 'GVECT_NUMBER', 'SMOOTH_GVECT_NUMBER']:
            parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)

        tagname = 'GAMMA_ONLY'
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

        tagname = 'SMALLBOX_FFT_GRID'
        fft_grid = []
        for attrname in ['nr1b', 'nr2b', 'nr3b']:
            fft_grid.append(parse_xml_child_attribute_int(tagname, attrname, target_tags))
        parsed_data[tagname.lower()] = fft_grid

    return parsed_data


def xml_card_symmetries(parsed_data, dom):
    cardname = 'SYMMETRIES'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['NUMBER_OF_SYMMETRIES', 'NUMBER_OF_BRAVAIS_SYMMETRIES']:
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_integer(tagname,target_tags)

    for tagname in ['INVERSION_SYMMETRY', 'DO_NOT_USE_TIME_REVERSAL', 'TIME_REVERSAL_FLAG', 'NO_TIME_REV_OPERATIONS']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    tagname = 'UNITS_FOR_SYMMETRIES'
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if metric not in ['crystal']:
        raise QEOutputParsingError(f'Error parsing attribute {attrname},' + \
                                   f' tag {tagname} inside ' + \
                                   f'{target_tags.tagName}, units unknown' )
    parsed_data['symmetries' + units_suffix] = metric

    # parse the symmetry matrices
    parsed_data['symmetries'] = []
    find_sym = True
    i = 0
    while find_sym:
        try:
            i += 1
            current_sym = {}
            tagname = 'SYMM.' + str(i)
            #a=target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            tagname2 = 'INFO'
            b = a.getElementsByTagName(tagname2)[0]
            attrname = 'NAME'
            value = str(b.getAttribute(attrname)).rstrip().replace('\n', '')
            current_sym['name'] = value

            try:
                attrname = 'T_REV'
                value = str(b.getAttribute(attrname)).rstrip().replace('\n', '')
                current_sym[attrname.lower()] = value
            except Exception:
                pass

            tagname2 = 'ROTATION'
            b = a.getElementsByTagName(tagname2)[0]
            c = [int(s) for s in b.childNodes[0].data.split()]
            current_sym[tagname2.lower()] = convert_list_to_matrix(c, 3, 3)

            for tagname2 in ['FRACTIONAL_TRANSLATION', 'EQUIVALENT_IONS']:  # not always present
                try:
                    b = a.getElementsByTagName(tagname2)[0]
                    if tagname2 == 'FRACTIONAL_TRANSLATION':
                        value = [float(s) for s in b.childNodes[0].data.split()]
                    else:
                        value = [int(s) for s in b.childNodes[0].data.split()]
                    current_sym[tagname2.lower()] = value
                except Exception:
                    raise

            parsed_data['symmetries'].append(current_sym)
        except IndexError:  # SYMM.i out of index
            find_sym = False

    return parsed_data


def xml_card_exchangecorrelation(parsed_data, dom):
    cardname = 'EXCHANGE_CORRELATION'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'DFT'
    parsed_data[(tagname+'_exchange_correlation').lower()] = \
        parse_xml_child_str(tagname,target_tags)

    tagname = 'LDA_PLUS_U_CALCULATION'
    try:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)
    except Exception:
        parsed_data[tagname.lower()] = False

    if parsed_data[tagname.lower()]:  # if it is a plus U calculation, I expect more infos
        tagname = 'HUBBARD_L'
        try:
            #a = target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            b = a.childNodes[0]
            c = b.data.replace('\n', '').split()
            value = [int(i) for i in c]
            parsed_data[tagname.lower()] = value
        except Exception:
            raise QEOutputParsingError('Error parsing tag '+\
                                       f'{tagname} inside {target_tags.tagName}.' )

        for tagname in ['HUBBARD_U', 'HUBBARD_ALPHA', 'HUBBARD_BETA', 'HUBBARD_J0']:
            try:
                #a = target_tags.getElementsByTagName(tagname)[0]
                a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
                b = a.childNodes[0]
                c = b.data.replace('\n', ' ').split()  # note the need of a white space!
                value = [float(i) * CONSTANTS.ry_to_ev for i in c]
                parsed_data[tagname.lower()] = value
            except Exception:
                raise QEOutputParsingError('Error parsing tag '+\
                                           f'{tagname} inside {target_tags.tagName}.')

        tagname = 'LDA_PLUS_U_KIND'
        try:
            parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)
        except Exception:
            pass

        tagname = 'U_PROJECTION_TYPE'
        try:
            parsed_data[tagname.lower()] = parse_xml_child_str(tagname, target_tags)
        except Exception:
            pass

        tagname = 'HUBBARD_J'
        try:
            #a=target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            b = a.childNodes[0]
            c = b.data.replace('\n', '').split()
            parsed_data[tagname.lower()] = convert_list_to_matrix(c, 3, 3)
        except Exception:
            pass

    try:
        tagname = 'NON_LOCAL_DF'
        parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)
    except Exception:
        pass

    try:
        tagname = 'VDW_KERNEL_NAME'
        parsed_data[tagname.lower()] = parse_xml_child_str(tagname, target_tags)
    except Exception:
        pass

    return parsed_data
