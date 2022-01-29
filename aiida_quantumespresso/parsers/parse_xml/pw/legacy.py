# -*- coding: utf-8 -*-
"""Code that was written to parse the legacy XML format of Quantum ESPRESSO, which was deprecated in version 6.4."""
import os
from xml.dom.minidom import parse, parseString

from qe_tools import CONSTANTS

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_xml.legacy import (
    parse_xml_child_attribute_str,
    parse_xml_child_bool,
    parse_xml_child_float,
    parse_xml_child_integer,
    read_xml_card,
    xml_card_cell,
    xml_card_exchangecorrelation,
    xml_card_header,
    xml_card_ions,
    xml_card_planewaves,
    xml_card_spin,
    xml_card_symmetries,
)
from aiida_quantumespresso.utils.mapping import get_logging_container

units_suffix = '_units'
default_energy_units = 'eV'
default_k_points_units = '1 / angstrom'
default_length_units = 'Angstrom'


def parse_pw_xml_pre_6_2(xml_file, dir_with_bands):
    """Parse the content of XML output file written by `pw.x` with the old schema-less XML format.

    :param xml_file: filelike object to the XML output file
    :param dir_with_bands: absolute filepath to directory containing k-point XML files
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    import copy
    from xml.parsers.expat import ExpatError

    logs = get_logging_container()

    # NOTE : I often assume that if the xml file has been written, it has no internal errors.
    try:
        dom = parse(xml_file)
    except ExpatError:
        logs.error.append('Error in XML parseString: bad format')
        parsed = {
            'bands': {},
            'structure': {},
        }
        return parsed, logs

    parsed_data = {}

    structure_dict = {}
    # CARD CELL
    structure_dict, lattice_vectors, volume = copy.deepcopy(xml_card_cell(structure_dict, dom))

    # CARD IONS
    structure_dict = copy.deepcopy(xml_card_ions(structure_dict, dom, lattice_vectors, volume))

    #CARD HEADER
    parsed_data = copy.deepcopy(xml_card_header(parsed_data, dom))

    # CARD CONTROL
    cardname = 'CONTROL'
    target_tags = read_xml_card(dom, cardname)
    for tagname in ['PP_CHECK_FLAG', 'LKPOINT_DIR', 'Q_REAL_SPACE', 'BETA_REAL_SPACE']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    # TODO: why this one isn't working? What is it actually?
#    # CARD MOVING_CELL
#
#    try:
#        target_tags = dom.getElementsByTagName('MOVING_CELL')[0]
#    except:
#        raise IOError
#
#    tagname='CELL_FACTOR'
#    parsed_data[tagname.lower()]=parse_xml_child_float(tagname,target_tags)

# CARD ELECTRIC_FIELD
    cardname = 'ELECTRIC_FIELD'
    target_tags = read_xml_card(dom, cardname)
    for tagname in ['HAS_ELECTRIC_FIELD', 'HAS_DIPOLE_CORRECTION']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)

    if parsed_data['has_electric_field'] or parsed_data['has_dipole_correction']:
        tagname = 'FIELD_DIRECTION'
        parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)

        for tagname in ['MAXIMUM_POSITION', 'INVERSE_REGION', 'FIELD_AMPLITUDE']:
            parsed_data[tagname.lower()] = parse_xml_child_float(tagname, target_tags)

    # CARD PLANE_WAVES
    parsed_data = copy.deepcopy(xml_card_planewaves(parsed_data, dom, 'pw'))

    # CARD SPIN
    parsed_data = copy.deepcopy(xml_card_spin(parsed_data, dom))

    # CARD BRILLOUIN ZONE
    cardname = 'BRILLOUIN_ZONE'
    target_tags = read_xml_card(dom, cardname)

    tagname = 'NUMBER_OF_K-POINTS'
    parsed_data[tagname.replace('-', '_').lower()] = parse_xml_child_integer(tagname, target_tags)

    tagname = 'UNITS_FOR_K-POINTS'
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if metric not in ['2 pi / a']:
        raise QEOutputParsingError(f'Error parsing attribute {attrname},' + \
                f' tag {tagname} inside {target_tags.tagName}, units unknown' )
    k_points_units = metric

    for tagname, param in [['MONKHORST_PACK_GRID', 'nk'], ['MONKHORST_PACK_OFFSET', 'k']]:
        try:
            #a = target_tags.getElementsByTagName(tagname)[0]
            a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            value = [int(a.getAttribute(param + str(i + 1))) for i in range(3)]
            parsed_data[tagname.replace('-', '_').lower()] = value
        except Exception:  # I might not use the monkhorst pack grid
            pass

    kpoints = []
    kpoints_weights = []

    tagname_prefix = 'K-POINT.'
    a_dict = {_.nodeName: _ for _ in target_tags.childNodes if _.nodeName.startswith(tagname_prefix)}

    try:
        import numpy
        for i in range(parsed_data['number_of_k_points']):
            tagname = f'{tagname_prefix}{i + 1}'
            #a = target_tags.getElementsByTagName(tagname)[0]
            a = a_dict[tagname]
            b = a.getAttribute('XYZ').replace('\n', '').rsplit()
            value = [float(s) for s in b]
            metric = k_points_units
            if metric == '2 pi / a':
                value = [2. * numpy.pi * float(s) / structure_dict['lattice_parameter'] for s in value]
                weight = float(a.getAttribute('WEIGHT'))
                kpoints.append(value)
                kpoints_weights.append(weight)
        parsed_data['k_points'] = kpoints
        parsed_data['k_points' + units_suffix] = default_k_points_units
        parsed_data['k_points_weights'] = kpoints_weights
    except Exception:
        raise QEOutputParsingError(f'Error parsing tag K-POINT.{i + 1} inside {target_tags.tagName}.')

    # I skip this card until someone will have a need for this.
#     try:
#         tagname='STARTING_K-POINTS'
#         num_starting_k_points=parse_xml_child_integer(tagname,target_tags)
#         # raise exception if there is no such a key
#         parsed_data[tagname.replace('-','_').lower()]=num_starting_k_points
#
#         if parsed_data.get('starting_k_points'):
#             try:
#                 kpoints=[]
#                 for i in range(parsed_data['starting_k_points']):
#                     tagname='K-POINT_START.'+str(i+1)
#                     a=target_tags.getElementsByTagName(tagname)[0]
#                     b=a.getAttribute('XYZ').replace('\n','').rsplit()
#                     value=[ float(s) for s in b ]
#                     metric=parsed_data['k_points_units']
#                     if metric=='2 pi / a':
#                         value=[ float(s)/parsed_data['lattice_parameter'] for s in value ]
#
#                         weight=float(a.getAttribute('WEIGHT'))
#
#                         kpoints.append([value,weight])
#
#                 parsed_data['k_point_start']=kpoints
#             except Exception:
#                 raise QEOutputParsingError('Error parsing tag {}'.format(tagname)+\
#                                            ' inside {}.'.format(target_tags.tagName ) )
#     except Exception:
#         if not parsed_data.get('starting_k_points'):
#             pass
#         else:
#             parsed_data['xml_warnings'].append("Warning: could not parse {}".format(tagname))

# tagname='NORM-OF-Q'
# TODO: decide if save this parameter
# parsed_data[tagname.replace('-','_').lower()]=parse_xml_child_float(tagname,target_tags)

# CARD BAND STRUCTURE INFO
    cardname = 'BAND_STRUCTURE_INFO'
    target_tags = read_xml_card(dom, cardname)

    for tagname in ['NUMBER_OF_SPIN_COMPONENTS', 'NUMBER_OF_ATOMIC_WFC', 'NUMBER_OF_BANDS']:
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_integer(tagname,target_tags)

    tagname = 'NON-COLINEAR_CALCULATION'
    parsed_data[tagname.replace('-','_').lower()] = \
        parse_xml_child_bool(tagname,target_tags)

    tagname = 'NUMBER_OF_ELECTRONS'
    parsed_data[tagname.replace('-','_').lower()] = \
        parse_xml_child_float(tagname,target_tags)

    tagname = 'UNITS_FOR_ENERGIES'
    attrname = 'UNITS'
    units = parse_xml_child_attribute_str(tagname, attrname, target_tags)
    if units not in ['hartree']:
        raise QEOutputParsingError(f"Expected energy units in Hartree. Got instead {parsed_data['energy_units']}")

    try:
        tagname = 'TWO_FERMI_ENERGIES'
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)
    except Exception:
        pass

    if parsed_data.get('two_fermi_energies', False):
        tagname = 'FERMI_ENERGY_UP'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * CONSTANTS.hartree_to_ev
        parsed_data[tagname.lower() + units_suffix] = default_energy_units
        tagname = 'FERMI_ENERGY_DOWN'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * CONSTANTS.hartree_to_ev
        parsed_data[tagname.lower() + units_suffix] = default_energy_units
    else:
        tagname = 'FERMI_ENERGY'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * CONSTANTS.hartree_to_ev
        parsed_data[tagname.lower() + units_suffix] = default_energy_units

    #CARD MAGNETIZATION_INIT
    cardname = 'MAGNETIZATION_INIT'
    target_tags = read_xml_card(dom, cardname)

    # 0 if false
    tagname = 'CONSTRAINT_MAG'
    parsed_data[tagname.lower()] = parse_xml_child_integer(tagname, target_tags)

    vec1 = []
    vec2 = []
    vec3 = []
    for i in range(structure_dict['number_of_species']):
        tagname = 'SPECIE.' + str(i + 1)
        #a=target_tags.getElementsByTagName(tagname)[0]
        a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        tagname2 = 'STARTING_MAGNETIZATION'
        vec1.append(parse_xml_child_float(tagname2, a))
        tagname2 = 'ANGLE1'
        vec2.append(parse_xml_child_float(tagname2, a))
        tagname2 = 'ANGLE2'
        vec3.append(parse_xml_child_float(tagname2, a))
    parsed_data['starting_magnetization'] = vec1
    parsed_data['magnetization_angle1'] = vec2
    parsed_data['magnetization_angle2'] = vec3

    #CARD OCCUPATIONS
    cardname = 'OCCUPATIONS'
    target_tags = read_xml_card(dom, cardname)
    for tagname in ['SMEARING_METHOD', 'TETRAHEDRON_METHOD', 'FIXED_OCCUPATIONS']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname, target_tags)
    if parsed_data['smearing_method']:
        parsed_data['occupations'] = 'smearing'
    elif parsed_data['tetrahedron_method']:
        parsed_data['occupations'] = 'tetrahedra'  # TODO: might also be tetrahedra_lin or tetrahedra_opt: check input?
    elif parsed_data['fixed_occupations']:
        parsed_data['occupations'] = 'fixed'

    # Remove the following deprecated keys
    for tagname in ['SMEARING_METHOD', 'TETRAHEDRON_METHOD', 'FIXED_OCCUPATIONS']:
        parsed_data.pop(tagname.lower())

    #CARD CHARGE-DENSITY
    cardname = 'CHARGE-DENSITY'
    target_tags = read_xml_card(dom, cardname)
    try:
        attrname = 'iotk_link'
        value = str(target_tags.getAttribute(attrname)).rstrip().replace('\n', '').lower()
        parsed_data[cardname.lower().rstrip().replace('-', '_')] = value
    except Exception:
        raise QEOutputParsingError(f'Error parsing attribute {attrname},' + \
                                   f' card {cardname}.')

    #CARD EIGENVALUES
    # Note: if this card is parsed, the dimension of the database grows very much!
    cardname = 'EIGENVALUES'
    target_tags = read_xml_card(dom, cardname)
    bands_dict = {}
    if dir_with_bands:
        try:
            occupations1 = []
            occupations2 = []
            bands1 = []
            bands2 = []
            for i in range(parsed_data['number_of_k_points']):
                tagname = 'K-POINT.' + str(i + 1)
                #a=target_tags.getElementsByTagName(tagname)[0]
                a = [_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

                def read_bands_and_occupations(eigenval_n):
                    # load the eigenval.xml file
                    with open(eigenval_n, 'r') as eigenval_f:
                        f = eigenval_f.read()

                    eig_dom = parseString(f)

                    tagname = 'UNITS_FOR_ENERGIES'
                    a = eig_dom.getElementsByTagName(tagname)[0]
                    attrname = 'UNITS'
                    metric = str(a.getAttribute(attrname))
                    if metric not in ['Hartree']:
                        raise QEOutputParsingError('Error parsing eigenvalues xml file, ' + \
                                                   f'units {metric} not implemented.')

                    tagname = 'EIGENVALUES'
                    a = eig_dom.getElementsByTagName(tagname)[0]
                    b = a.childNodes[0]
                    value_e = [float(s) * CONSTANTS.hartree_to_ev for s in b.data.split()]

                    tagname = 'OCCUPATIONS'
                    a = eig_dom.getElementsByTagName(tagname)[0]
                    b = a.childNodes[0]
                    value_o = [float(s) for s in b.data.split()]
                    return value_e, value_o

                # two cases: in cases of magnetic calculations, I have both spins
                try:
                    tagname2 = 'DATAFILE'
                    b = a.getElementsByTagName(tagname2)[0]
                    attrname = 'iotk_link'
                    value = str(b.getAttribute(attrname)).rstrip().replace('\n', '')
                    eigenval_n = os.path.join(dir_with_bands, value)

                    value_e, value_o = read_bands_and_occupations(eigenval_n)
                    bands1.append(value_e)
                    occupations1.append(value_o)

                except IndexError:
                    tagname2 = 'DATAFILE.1'
                    b1 = a.getElementsByTagName(tagname2)[0]
                    tagname2 = 'DATAFILE.2'
                    b2 = a.getElementsByTagName(tagname2)[0]
                    attrname = 'iotk_link'
                    value1 = str(b1.getAttribute(attrname)).rstrip().replace('\n', '')
                    value2 = str(b2.getAttribute(attrname)).rstrip().replace('\n', '')

                    eigenval_n = os.path.join(dir_with_bands, value1)
                    value_e, value_o = read_bands_and_occupations(eigenval_n)
                    bands1.append(value_e)
                    occupations1.append(value_o)

                    eigenval_n = os.path.join(dir_with_bands, value2)
                    value_e, value_o = read_bands_and_occupations(eigenval_n)
                    bands2.append(value_e)
                    occupations2.append(value_o)

            occupations = [occupations1]
            bands = [bands1]
            if occupations2:
                occupations.append(occupations2)
            if bands2:
                bands.append(bands2)

            bands_dict['occupations'] = occupations
            bands_dict['bands'] = bands
            bands_dict['bands' + units_suffix] = default_energy_units
        except Exception as exception:
            raise QEOutputParsingError(f'Error parsing card {tagname}: {exception.__class__.__name__} {exception}')


#     if dir_with_bands:
#         # if there is at least an empty band:
#         if parsed_data['smearing_method'] or  \
#            parsed_data['number_of_electrons']/2. < parsed_data['number_of_bands']:
#
#             #TODO: currently I do it only for non magnetic systems
#             if len(bands_dict['occupations'])==1:
#             # initialize lumo
#                 lumo = parsed_data['homo']+10000.0
#                 for list_bands in bands_dict['bands']:
#                     for value in list_bands:
#                         if (value > parsed_data['fermi_energy']) and (value<lumo):
#                             lumo=value
#                 if (lumo==parsed_data['homo']+10000.0) or lumo<=parsed_data['fermi_energy']:
#                     #might be an error for bandgap larger than 10000 eV...
#                     raise QEOutputParsingError('Error while searching for LUMO.')
#                 parsed_data['lumo']=lumo
#                 parsed_data['lumo'+units_suffix] = default_energy_units

# CARD symmetries
    parsed_data = copy.deepcopy(xml_card_symmetries(parsed_data, dom))

    # CARD EXCHANGE_CORRELATION
    parsed_data = copy.deepcopy(xml_card_exchangecorrelation(parsed_data, dom))

    parsed_data['bands'] = bands_dict
    parsed_data['structure'] = structure_dict

    return parsed_data, logs
