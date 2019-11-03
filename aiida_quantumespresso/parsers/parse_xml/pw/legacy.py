# -*- coding: utf-8 -*-
"""Code that was written to parse the legacy XML format of Quantum ESPRESSO, which was deprecated in version 6.4."""
from __future__ import absolute_import
from __future__ import print_function

import os
import string
from xml.dom.minidom import parse, parseString, Element
from six.moves import range

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.utils.mapping import get_logging_container
from qe_tools.constants import ry_to_ev, hartree_to_ev, bohr_to_ang

units_suffix = '_units'
default_energy_units = 'eV'
default_k_points_units = '1 / angstrom'
default_length_units = 'Angstrom'


def parse_pw_xml_pre_6_2(xml_file, dir_with_bands, include_deprecated_keys=False):
    """Parse the content of XML output file written by `pw.x` with the old schema-less XML format.

    :param xml_file: filelike object to the XML output file
    :param dir_with_bands: absolute filepath to directory containing k-point XML files
    :param include_deprecated_v2_keys: boolean, if True, includes deprecated keys from old parser v2
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
    structure_dict,lattice_vectors,volume = copy.deepcopy(xml_card_cell(structure_dict,dom))

    # CARD IONS
    structure_dict = copy.deepcopy(xml_card_ions(structure_dict,dom,lattice_vectors,volume))

    #CARD HEADER
    parsed_data = copy.deepcopy(xml_card_header(parsed_data, dom))

    # CARD CONTROL
    cardname='CONTROL'
    target_tags=read_xml_card(dom,cardname)
    for tagname in ['PP_CHECK_FLAG','LKPOINT_DIR',
                    'Q_REAL_SPACE','BETA_REAL_SPACE']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname,target_tags)

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
    cardname='ELECTRIC_FIELD'
    target_tags=read_xml_card(dom,cardname)
    for tagname in ['HAS_ELECTRIC_FIELD','HAS_DIPOLE_CORRECTION']:
        parsed_data[tagname.lower()]=parse_xml_child_bool(tagname,target_tags)

    if parsed_data['has_electric_field'] or parsed_data['has_dipole_correction']:
        tagname='FIELD_DIRECTION'
        parsed_data[tagname.lower()]=parse_xml_child_integer(tagname,target_tags)

        for tagname in ['MAXIMUM_POSITION','INVERSE_REGION','FIELD_AMPLITUDE']:
            parsed_data[tagname.lower()]=parse_xml_child_float(tagname,target_tags)

    # CARD PLANE_WAVES
    parsed_data = copy.deepcopy(xml_card_planewaves(parsed_data,dom,'pw'))

    # CARD SPIN
    parsed_data = copy.deepcopy(xml_card_spin(parsed_data,dom))

    # CARD BRILLOUIN ZONE
    cardname='BRILLOUIN_ZONE'
    target_tags=read_xml_card(dom,cardname)

    tagname='NUMBER_OF_K-POINTS'
    parsed_data[tagname.replace('-','_').lower()] = parse_xml_child_integer(tagname,target_tags)

    tagname = 'UNITS_FOR_K-POINTS'
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname,attrname,target_tags)
    if metric not in ['2 pi / a']:
        raise QEOutputParsingError('Error parsing attribute {},'.format(attrname) + \
                ' tag {} inside {}, units unknown'.format(tagname, target_tags.tagName) )
    k_points_units = metric

    for tagname, param in [ ['MONKHORST_PACK_GRID','nk'],['MONKHORST_PACK_OFFSET','k'] ]:
        try:
            #a = target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            value = [int(a.getAttribute(param+str(i+1))) for i in range(3)]
            parsed_data[tagname.replace('-','_').lower()] = value
        except Exception: # I might not use the monkhorst pack grid
            pass


    kpoints = []
    kpoints_weights = []

    tagname_prefix = 'K-POINT.'
    a_dict={_.nodeName: _ for _ in target_tags.childNodes
            if _.nodeName.startswith(tagname_prefix)}

    try:
        import numpy
        for i in range(parsed_data['number_of_k_points']):
            tagname = '{}{}'.format(tagname_prefix,i+1)
            #a = target_tags.getElementsByTagName(tagname)[0]
            a=a_dict[tagname]
            b = a.getAttribute('XYZ').replace('\n','').rsplit()
            value = [ float(s) for s in b ]
            metric = k_points_units
            if metric=='2 pi / a':
                value = [ 2.*numpy.pi*float(s)/structure_dict['lattice_parameter'] for s in value ]
                weight = float(a.getAttribute('WEIGHT'))
                kpoints.append(value)
                kpoints_weights.append(weight)
        parsed_data['k_points']=kpoints
        parsed_data['k_points'+units_suffix] = default_k_points_units
        parsed_data['k_points_weights'] = kpoints_weights
    except Exception:
        raise QEOutputParsingError('Error parsing tag K-POINT.{} inside {}.'.format(i+1,target_tags.tagName) )


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
    target_tags = read_xml_card(dom,cardname)

    for tagname in ['NUMBER_OF_SPIN_COMPONENTS','NUMBER_OF_ATOMIC_WFC','NUMBER_OF_BANDS']:
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_integer(tagname,target_tags)

    tagname='NON-COLINEAR_CALCULATION'
    parsed_data[tagname.replace('-','_').lower()] = \
        parse_xml_child_bool(tagname,target_tags)

    tagname='NUMBER_OF_ELECTRONS'
    parsed_data[tagname.replace('-','_').lower()] = \
        parse_xml_child_float(tagname,target_tags)

    tagname = 'UNITS_FOR_ENERGIES'
    attrname = 'UNITS'
    units = parse_xml_child_attribute_str(tagname,attrname,target_tags)
    if units not in ['hartree']:
        raise QEOutputParsingError('Expected energy units in Hartree. Got instead {}'.format(parsed_data['energy_units']))

    try:
        tagname='TWO_FERMI_ENERGIES'
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname,target_tags)
    except Exception:
        pass

    if parsed_data.get('two_fermi_energies',False):
        tagname = 'FERMI_ENERGY_UP'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * hartree_to_ev
        parsed_data[tagname.lower()+units_suffix] = default_energy_units
        tagname = 'FERMI_ENERGY_DOWN'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * hartree_to_ev
        parsed_data[tagname.lower()+units_suffix] = default_energy_units
    else:
        tagname = 'FERMI_ENERGY'
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_float(tagname,target_tags) * hartree_to_ev
        parsed_data[tagname.lower()+units_suffix] = default_energy_units

    #CARD MAGNETIZATION_INIT
    cardname = 'MAGNETIZATION_INIT'
    target_tags = read_xml_card(dom,cardname)

    # 0 if false
    tagname='CONSTRAINT_MAG'
    parsed_data[tagname.lower()] = parse_xml_child_integer(tagname,target_tags)

    vec1 = []
    vec2 = []
    vec3 = []
    for i in range(structure_dict['number_of_species']):
        tagname='SPECIE.'+str(i+1)
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        tagname2='STARTING_MAGNETIZATION'
        vec1.append(parse_xml_child_float(tagname2,a))
        tagname2='ANGLE1'
        vec2.append(parse_xml_child_float(tagname2,a))
        tagname2='ANGLE2'
        vec3.append(parse_xml_child_float(tagname2,a))
    parsed_data['starting_magnetization'] = vec1
    parsed_data['magnetization_angle1'] = vec2
    parsed_data['magnetization_angle2'] = vec3

    #CARD OCCUPATIONS
    cardname = 'OCCUPATIONS'
    target_tags = read_xml_card(dom,cardname)
    for tagname in ['SMEARING_METHOD','TETRAHEDRON_METHOD','FIXED_OCCUPATIONS']:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname,target_tags)
    if parsed_data['smearing_method']:
        parsed_data['occupations'] = 'smearing'
    elif parsed_data['tetrahedron_method']:
        parsed_data['occupations'] = 'tetrahedra'  # TODO: might also be tetrahedra_lin or tetrahedra_opt: check input?
    elif parsed_data['fixed_occupations']:
        parsed_data['occupations'] = 'fixed'
    if not include_deprecated_keys:
        for tagname in ['SMEARING_METHOD','TETRAHEDRON_METHOD','FIXED_OCCUPATIONS']:
            parsed_data.pop(tagname.lower())

    #CARD CHARGE-DENSITY
    cardname='CHARGE-DENSITY'
    target_tags=read_xml_card(dom,cardname)
    try:
        attrname='iotk_link'
        value=str(target_tags.getAttribute(attrname)).rstrip().replace('\n','').lower()
        parsed_data[cardname.lower().rstrip().replace('-','_')]=value
    except Exception:
        raise QEOutputParsingError('Error parsing attribute {},'.format(attrname) + \
                                   ' card {}.'.format(cardname))

    #CARD EIGENVALUES
    # Note: if this card is parsed, the dimension of the database grows very much!
    cardname='EIGENVALUES'
    target_tags=read_xml_card(dom,cardname)
    bands_dict = {}
    if dir_with_bands:
        try:
            occupations1 = []
            occupations2 = []
            bands1 = []
            bands2 = []
            for i in range(parsed_data['number_of_k_points']):
                tagname='K-POINT.'+str(i+1)
                #a=target_tags.getElementsByTagName(tagname)[0]
                a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

                def read_bands_and_occupations(eigenval_n):
                    # load the eigenval.xml file
                    with open(eigenval_n,'r') as eigenval_f:
                        f = eigenval_f.read()

                    eig_dom = parseString(f)

                    tagname = 'UNITS_FOR_ENERGIES'
                    a = eig_dom.getElementsByTagName(tagname)[0]
                    attrname = 'UNITS'
                    metric = str(a.getAttribute(attrname))
                    if metric not in ['Hartree']:
                        raise QEOutputParsingError('Error parsing eigenvalues xml file, ' + \
                                                   'units {} not implemented.'.format(metric))

                    tagname='EIGENVALUES'
                    a=eig_dom.getElementsByTagName(tagname)[0]
                    b=a.childNodes[0]
                    value_e = [ float(s)*hartree_to_ev for s in b.data.split() ]

                    tagname='OCCUPATIONS'
                    a = eig_dom.getElementsByTagName(tagname)[0]
                    b = a.childNodes[0]
                    value_o = [ float(s) for s in b.data.split() ]
                    return value_e,value_o

                # two cases: in cases of magnetic calculations, I have both spins
                try:
                    tagname2 = 'DATAFILE'
                    b = a.getElementsByTagName(tagname2)[0]
                    attrname = 'iotk_link'
                    value = str(b.getAttribute(attrname)).rstrip().replace('\n','')
                    eigenval_n =  os.path.join(dir_with_bands,value)

                    value_e,value_o = read_bands_and_occupations(eigenval_n)
                    bands1.append(value_e)
                    occupations1.append(value_o)

                except IndexError:
                    tagname2='DATAFILE.1'
                    b1 = a.getElementsByTagName(tagname2)[0]
                    tagname2='DATAFILE.2'
                    b2 = a.getElementsByTagName(tagname2)[0]
                    attrname = 'iotk_link'
                    value1 = str(b1.getAttribute(attrname)).rstrip().replace('\n','')
                    value2 = str(b2.getAttribute(attrname)).rstrip().replace('\n','')

                    eigenval_n =  os.path.join(dir_with_bands,value1)
                    value_e,value_o = read_bands_and_occupations(eigenval_n)
                    bands1.append(value_e)
                    occupations1.append(value_o)

                    eigenval_n =  os.path.join(dir_with_bands,value2)
                    value_e,value_o = read_bands_and_occupations(eigenval_n)
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
            bands_dict['bands'+units_suffix] = default_energy_units
        except Exception as exception:
            raise QEOutputParsingError('Error parsing card {}: {} {}'.format(
                tagname,exception.__class__.__name__,exception
                ))

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
    parsed_data = copy.deepcopy(xml_card_symmetries(parsed_data,dom))

    # CARD EXCHANGE_CORRELATION
    parsed_data = copy.deepcopy(xml_card_exchangecorrelation(parsed_data,dom))

    parsed_data['bands'] = bands_dict
    parsed_data['structure'] = structure_dict

    return parsed_data, logs


def cell_volume(a1,a2,a3):
    r"""
    returns the volume of the primitive cell: :math:`|\vec a_1\cdot(\vec a_2\cross \vec a_3)|`
    """
    a_mid_0 = a2[1]*a3[2] - a2[2]*a3[1]
    a_mid_1 = a2[2]*a3[0] - a2[0]*a3[2]
    a_mid_2 = a2[0]*a3[1] - a2[1]*a3[0]
    return abs(float(a1[0]*a_mid_0 + a1[1]*a_mid_1 + a1[2]*a_mid_2))


# In the following, some functions that helps the parsing of
# the xml file of QE v5.0.x (version below not tested)
def read_xml_card(dom,cardname):
    try:
        root_node = [_ for _ in dom.childNodes if
                    isinstance(_, Element)
                        and _.nodeName == 'Root'][0]
        the_card = [_ for _ in root_node.childNodes if _.nodeName == cardname][0]
        #the_card = dom.getElementsByTagName(cardname)[0]
        return the_card
    except Exception as e:
        print(e)
        raise QEOutputParsingError('Error parsing tag {}'.format(cardname) )

def parse_xml_child_integer(tagname,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.childNodes[0]
        return int(b.data)
    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}'
                                   .format(tagname,target_tags.tagName) )

def parse_xml_child_float(tagname,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.childNodes[0]
        return float(b.data)
    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}'\
                                   .format(tagname, target_tags.tagName ) )

def parse_xml_child_bool(tagname,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.childNodes[0]
        return str2bool(b.data)
    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}'\
                                   .format(tagname, target_tags.tagName) )

def str2bool(string):
    try:
        false_items=['f','0','false','no']
        true_items=['t','1','true','yes']
        string=str(string.lower().strip())
        if string in false_items:
            return False
        if string in true_items:
            return True
        else:
            raise QEOutputParsingError('Error converting string '
                                       '{} to boolean value.'.format(string) )
    except Exception:
        raise QEOutputParsingError('Error converting string to boolean.')

def parse_xml_child_str(tagname,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.childNodes[0]
        return str(b.data).rstrip().replace('\n','')
    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}'\
                                   .format(tagname, target_tags.tagName) )

def parse_xml_child_attribute_str(tagname,attributename,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        value=str(a.getAttribute(attributename))
        return value.rstrip().replace('\n','').lower()
    except Exception:
        raise QEOutputParsingError('Error parsing attribute {}, tag {} inside {}'
                                   .format(attributename,tagname,target_tags.tagName) )

def parse_xml_child_attribute_int(tagname,attributename,target_tags):
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        value=int(a.getAttribute(attributename))
        return value
    except Exception:
        raise QEOutputParsingError('Error parsing attribute {}, tag {} inside {}'
                                    .format(attributename,tagname,target_tags.tagName) )


def convert_list_to_matrix(in_matrix,n_rows,n_columns):
    """converts a list into a list of lists (a matrix like) with n_rows and n_columns."""
    return [ in_matrix[j:j+n_rows] for j in range(0,n_rows*n_columns,n_rows) ]

def xml_card_cell(parsed_data,dom):
    #CARD CELL of QE output

    cardname = 'CELL'
    target_tags = read_xml_card(dom,cardname)

    for tagname in ['NON-PERIODIC_CELL_CORRECTION','BRAVAIS_LATTICE']:
        parsed_data[tagname.replace('-','_').lower()] = parse_xml_child_str(tagname,target_tags)

    tagname = 'LATTICE_PARAMETER'
    value = parse_xml_child_float(tagname,target_tags)
    parsed_data[tagname.replace('-','_').lower()+'_xml'] = value
    attrname = 'UNITS'
    metric = parse_xml_child_attribute_str(tagname,attrname,target_tags)
    if metric not in ['bohr','angstrom']:
        raise QEOutputParsingError('Error parsing attribute {}, tag {} inside {}, units not found'
                                   .format(attrname,tagname,target_tags.tagName) )
    if metric == 'bohr':
        value *= bohr_to_ang
    parsed_data[tagname.replace('-','_').lower()] = value

    tagname='CELL_DIMENSIONS'
    try:
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.childNodes[0]
        c=b.data.replace('\n','').split()
        value=[ float(i) for i in c ]
        parsed_data[tagname.replace('-','_').lower()]=value
    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}.'.format(tagname,target_tags.tagName) )

    tagname = 'DIRECT_LATTICE_VECTORS'
    lattice_vectors = []
    try:
        second_tagname='UNITS_FOR_DIRECT_LATTICE_VECTORS'
        #a=target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
        b=a.getElementsByTagName('UNITS_FOR_DIRECT_LATTICE_VECTORS')[0]
        value=str(b.getAttribute('UNITS')).lower()
        parsed_data[second_tagname.replace('-','_').lower()]=value

        metric = value
        if metric not in ['bohr','angstroms']: # REMEMBER TO CHECK THE UNITS AT THE END OF THE FUNCTION
            raise QEOutputParsingError('Error parsing tag {} inside {}: units not supported: {}'
                                       .format(tagname,target_tags.tagName,metric) )

        lattice_vectors = []
        for second_tagname in ['a1','a2','a3']:
            #b = a.getElementsByTagName(second_tagname)[0]
            b=[_ for _ in a.childNodes if _.nodeName == second_tagname][0]
            c = b.childNodes[0]
            d = c.data.replace('\n','').split()
            value = [ float(i) for i in d ]
            if metric=='bohr':
                value = [ bohr_to_ang*float(s) for s in value ]
            lattice_vectors.append(value)

        volume = cell_volume(lattice_vectors[0],lattice_vectors[1],lattice_vectors[2])

    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {} inside {}.'
                                   .format(tagname,target_tags.tagName,cardname) )
    # NOTE: lattice_vectors will be saved later together with card IONS.atom

    tagname = 'RECIPROCAL_LATTICE_VECTORS'
    try:
        #a = target_tags.getElementsByTagName(tagname)[0]
        a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

        second_tagname = 'UNITS_FOR_RECIPROCAL_LATTICE_VECTORS'
        b = a.getElementsByTagName(second_tagname)[0]
        value = str(b.getAttribute('UNITS')).lower()
        parsed_data[second_tagname.replace('-','_').lower()]=value

        metric=value
        # NOTE: output is given in 2 pi / a [ang ^ -1]
        if metric not in ['2 pi / a']:
            raise QEOutputParsingError('Error parsing tag {} inside {}: units {} not supported'
                                       .format(tagname,target_tags.tagName,metric) )

        # reciprocal_lattice_vectors
        this_matrix = []
        for second_tagname in ['b1','b2','b3']:
            b = a.getElementsByTagName(second_tagname)[0]
            c = b.childNodes[0]
            d = c.data.replace('\n','').split()
            value = [ float(i) for i in d ]
            if metric == '2 pi / a':
                value=[ float(s)/parsed_data['lattice_parameter'] for s in value ]
            this_matrix.append(value)
        parsed_data['reciprocal_lattice_vectors'] = this_matrix

    except Exception:
        raise QEOutputParsingError('Error parsing tag {} inside {}.'
                                   .format(tagname,target_tags.tagName) )
    return parsed_data,lattice_vectors,volume

def xml_card_ions(parsed_data,dom,lattice_vectors,volume):
    cardname='IONS'
    target_tags=read_xml_card(dom,cardname)

    for tagname in ['NUMBER_OF_ATOMS','NUMBER_OF_SPECIES']:
        parsed_data[tagname.lower()]=parse_xml_child_integer(tagname,target_tags)

    tagname='UNITS_FOR_ATOMIC_MASSES'
    attrname='UNITS'
    parsed_data[tagname.lower()]=parse_xml_child_attribute_str(tagname,attrname,target_tags)

    try:
        parsed_data['species']={}
        parsed_data['species']['index'] =[]
        parsed_data['species']['type']  =[]
        parsed_data['species']['mass']  =[]
        parsed_data['species']['pseudo']=[]
        for i in range(parsed_data['number_of_species']):
            tagname='SPECIE.'+str(i+1)
            parsed_data['species']['index'].append(i+1)

            #a=target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]

            tagname2='ATOM_TYPE'
            parsed_data['species']['type'].append(parse_xml_child_str(tagname2,a))

            tagname2='MASS'
            parsed_data['species']['mass'].append(parse_xml_child_float(tagname2,a))

            tagname2='PSEUDO'
            parsed_data['species']['pseudo'].append(parse_xml_child_str(tagname2,a))

        tagname='UNITS_FOR_ATOMIC_POSITIONS'
        attrname='UNITS'
        parsed_data[tagname.lower()]=parse_xml_child_attribute_str(tagname,attrname,target_tags)
    except:
        raise QEOutputParsingError('Error parsing tag SPECIE.# inside %s.'% (target_tags.tagName ) )
# TODO convert the units
# if parsed_data['units_for_atomic_positions'] not in ['alat','bohr','angstrom']:

    try:
        atomlist=[]
        atoms_index_list=[]
        atoms_if_pos_list=[]
        tagslist=[]
        for i in range(parsed_data['number_of_atoms']):
            tagname='ATOM.'+str(i+1)
            # USELESS AT THE MOMENT, I DON'T SAVE IT
            # parsed_data['atoms']['list_index']=i
            #a=target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            tagname2='INDEX'
            b=int(a.getAttribute(tagname2))
            atoms_index_list.append(b)
            tagname2='SPECIES'

            chem_symbol=str(a.getAttribute(tagname2)).rstrip().replace('\n','')
            # I check if it is a subspecie
            chem_symbol_digits = ''.join([i for i in chem_symbol if i in string.digits])
            try:
                tagslist.append(int(chem_symbol_digits))
            except ValueError:
                # If I can't parse the digit, it is probably not there: I add a None to the tagslist
                tagslist.append(None)
            # I remove the symbols
            chem_symbol = ''.join(i for i in chem_symbol if not i.isdigit())

            tagname2='tau'
            b = a.getAttribute(tagname2)
            tau = [float(s) for s in b.rstrip().replace('\n','').split()]
            metric = parsed_data['units_for_atomic_positions']
            if metric not in ['alat','bohr','angstrom']: # REMEMBER TO CONVERT AT THE END
                raise QEOutputParsingError('Error parsing tag %s inside %s'% (tagname, target_tags.tagName ) )
            if metric=='alat':
                tau=[ parsed_data['lattice_parameter_xml']*float(s) for s in tau ]
            elif metric=='bohr':
                tau=[ bohr_to_ang*float(s) for s in tau ]
            atomlist.append([chem_symbol,tau])
            tagname2='if_pos'
            b=a.getAttribute(tagname2)
            if_pos=[int(s) for s in b.rstrip().replace('\n','').split()]
            atoms_if_pos_list.append(if_pos)
        parsed_data['atoms']=atomlist
        parsed_data['atoms_index_list']=atoms_index_list
        parsed_data['atoms_if_pos_list']=atoms_if_pos_list
        cell={}
        cell['lattice_vectors']=lattice_vectors
        cell['volume']=volume
        cell['atoms']=atomlist
        cell['tagslist'] = tagslist
        parsed_data['cell']=cell
    except Exception:
        raise QEOutputParsingError('Error parsing tag ATOM.# inside %s.'% (target_tags.tagName ) )
    # saving data together with cell parameters. Did so for better compatibility with ASE.

    # correct some units that have been converted in
    parsed_data['atomic_positions'+units_suffix] = default_length_units
    parsed_data['direct_lattice_vectors'+units_suffix] = default_length_units

    return parsed_data

def xml_card_spin(parsed_data,dom):
    cardname='SPIN'
    target_tags=read_xml_card(dom,cardname)

    for tagname in ['LSDA','NON-COLINEAR_CALCULATION',
                    'SPIN-ORBIT_CALCULATION','SPIN-ORBIT_DOMAG']:
        parsed_data[tagname.replace('-','_').lower()
                    ] = parse_xml_child_bool(tagname,target_tags)

    return parsed_data

def xml_card_header(parsed_data,dom):
    cardname='HEADER'
    target_tags=read_xml_card(dom,cardname)

    for tagname in ['FORMAT','CREATOR']:
        for attrname in ['NAME','VERSION']:
            parsed_data[(tagname+'_'+attrname).lower()
                        ] = parse_xml_child_attribute_str(tagname,attrname,target_tags)

    return parsed_data

def xml_card_planewaves(parsed_data,dom,calctype):
    if calctype not in ['pw','cp']:
        raise ValueError("Input flag not accepted, must be 'cp' or 'pw'")

    cardname='PLANE_WAVES'
    target_tags=read_xml_card(dom,cardname)

    tagname = 'UNITS_FOR_CUTOFF'
    attrname = 'UNITS'
    units = parse_xml_child_attribute_str(tagname,attrname,target_tags).lower()
    if 'hartree' not in units:
        if 'rydberg' not in units:
            raise QEOutputParsingError('Units {} are not supported by parser'.format(units))
    else:
        if 'hartree' in units:
            conv_fac = hartree_to_ev
        else:
            conv_fac = ry_to_ev

        tagname='WFC_CUTOFF'
        parsed_data[tagname.lower()] = parse_xml_child_float(tagname,target_tags)*conv_fac
        parsed_data[tagname.lower()+units_suffix] = default_energy_units

        tagname='RHO_CUTOFF'
        parsed_data[tagname.lower()] = parse_xml_child_float(tagname,target_tags)*conv_fac
        parsed_data[tagname.lower()+units_suffix] = default_energy_units

    for tagname in [ 'FFT_GRID','SMOOTH_FFT_GRID' ]:
        grid = []
        for attrname in ['nr1','nr2','nr3']:
            if 'SMOOTH' in tagname:
                attrname += 's'
            grid.append(parse_xml_child_attribute_int(tagname,attrname,target_tags))
        parsed_data[tagname.lower()] = grid

    if calctype == 'cp':

        for tagname in ['MAX_NUMBER_OF_GK-VECTORS','GVECT_NUMBER','SMOOTH_GVECT_NUMBER' ]:
            parsed_data[tagname.lower()] = parse_xml_child_integer(tagname,target_tags)

        tagname='GAMMA_ONLY'
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname,target_tags)

        tagname='SMALLBOX_FFT_GRID'
        fft_grid = []
        for attrname in ['nr1b','nr2b','nr3b']:
            fft_grid.append(parse_xml_child_attribute_int(tagname,attrname,target_tags))
        parsed_data[tagname.lower()] = fft_grid

    return parsed_data

def xml_card_symmetries(parsed_data,dom):
    cardname='SYMMETRIES'
    target_tags=read_xml_card(dom,cardname)

    for tagname in ['NUMBER_OF_SYMMETRIES','NUMBER_OF_BRAVAIS_SYMMETRIES']:
        parsed_data[tagname.replace('-','_').lower()] = \
            parse_xml_child_integer(tagname,target_tags)

    for tagname in ['INVERSION_SYMMETRY','DO_NOT_USE_TIME_REVERSAL',
                    'TIME_REVERSAL_FLAG','NO_TIME_REV_OPERATIONS']:
        parsed_data[tagname.lower()]=parse_xml_child_bool(tagname,target_tags)

    tagname='UNITS_FOR_SYMMETRIES'
    attrname='UNITS'
    metric=parse_xml_child_attribute_str(tagname,attrname,target_tags)
    if metric not in ['crystal']:
        raise QEOutputParsingError('Error parsing attribute {},'.format(attrname) + \
                                   ' tag {} inside '.format(tagname) + \
                                   '{}, units unknown'.format(target_tags.tagName ) )
    parsed_data['symmetries'+units_suffix] = metric

    # parse the symmetry matrices
    parsed_data['symmetries']=[]
    find_sym=True
    i=0
    while find_sym:
        try:
            i+=1
            current_sym={}
            tagname='SYMM.'+str(i)
            #a=target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            tagname2='INFO'
            b=a.getElementsByTagName(tagname2)[0]
            attrname='NAME'
            value=str(b.getAttribute(attrname)).rstrip().replace('\n','')
            current_sym['name']=value

            try:
                attrname='T_REV'
                value=str(b.getAttribute(attrname)).rstrip().replace('\n','')
                current_sym[attrname.lower()]=value
            except Exception:
                pass

            tagname2='ROTATION'
            b=a.getElementsByTagName(tagname2)[0]
            c=[ int(s) for s in b.childNodes[0].data.split() ]
            current_sym[tagname2.lower()] = convert_list_to_matrix(c,3,3)

            for tagname2 in ['FRACTIONAL_TRANSLATION','EQUIVALENT_IONS']: # not always present
                try:
                    b = a.getElementsByTagName(tagname2)[0]
                    if tagname2 == 'FRACTIONAL_TRANSLATION':
                        value = [ float(s) for s in b.childNodes[0].data.split() ]
                    else:
                        value = [ int(s) for s in b.childNodes[0].data.split() ]
                    current_sym[tagname2.lower()] = value
                except Exception:
                    raise

            parsed_data['symmetries'].append(current_sym)
        except IndexError: # SYMM.i out of index
            find_sym=False

    return parsed_data

def xml_card_exchangecorrelation(parsed_data,dom):
    cardname='EXCHANGE_CORRELATION'
    target_tags=read_xml_card(dom,cardname)

    tagname='DFT'
    parsed_data[(tagname+'_exchange_correlation').lower()] = \
        parse_xml_child_str(tagname,target_tags)

    tagname='LDA_PLUS_U_CALCULATION'
    try:
        parsed_data[tagname.lower()] = parse_xml_child_bool(tagname,target_tags)
    except Exception:
        parsed_data[tagname.lower()] = False

    if parsed_data[tagname.lower()]: # if it is a plus U calculation, I expect more infos
        tagname = 'HUBBARD_L'
        try:
            #a = target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            b = a.childNodes[0]
            c = b.data.replace('\n','').split()
            value = [int(i) for i in c]
            parsed_data[tagname.lower()] = value
        except Exception:
            raise QEOutputParsingError('Error parsing tag '+\
                                       '{} inside {}.'.format(tagname, target_tags.tagName) )

        for tagname in ['HUBBARD_U','HUBBARD_ALPHA','HUBBARD_BETA','HUBBARD_J0']:
            try:
                #a = target_tags.getElementsByTagName(tagname)[0]
                a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
                b = a.childNodes[0]
                c = b.data.replace('\n',' ').split() # note the need of a white space!
                value = [float(i)*ry_to_ev for i in c]
                parsed_data[tagname.lower()] = value
            except Exception:
                raise QEOutputParsingError('Error parsing tag '+\
                                           '{} inside {}.'.format(tagname, target_tags.tagName))

        tagname = 'LDA_PLUS_U_KIND'
        try:
            parsed_data[tagname.lower()] = parse_xml_child_integer(tagname,target_tags)
        except Exception:
            pass

        tagname = 'U_PROJECTION_TYPE'
        try:
            parsed_data[tagname.lower()] = parse_xml_child_str(tagname,target_tags)
        except Exception:
            pass

        tagname = 'HUBBARD_J'
        try:
            #a=target_tags.getElementsByTagName(tagname)[0]
            a=[_ for _ in target_tags.childNodes if _.nodeName == tagname][0]
            b=a.childNodes[0]
            c=b.data.replace('\n','').split()
            parsed_data[tagname.lower()] = convert_list_to_matrix(c,3,3)
        except Exception:
            pass

    try:
        tagname='NON_LOCAL_DF'
        parsed_data[tagname.lower()] = parse_xml_child_integer(tagname,target_tags)
    except Exception:
        pass

    try:
        tagname='VDW_KERNEL_NAME'
        parsed_data[tagname.lower()] = parse_xml_child_str(tagname,target_tags)
    except Exception:
        pass

    return parsed_data
