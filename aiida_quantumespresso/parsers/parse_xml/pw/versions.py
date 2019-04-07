# -*- coding: utf-8 -*-
import enum
import os
from defusedxml import ElementTree

DEFAULT_SCHEMA_FILENAME = 'qes-1.0.xsd'


class QeXmlVersion(enum.Enum):
    """An enum with the versions of XML output file known to exist for Quantum ESPRESSO"""

    PRE_6_2 = 0
    POST_6_2 = 1


def get_xml_file_version(xml_file):
    """
    Get the version of the Quantum ESPRESSO pw.x XML output file

    :param xml_file: absolute path to the XML file
    :raises ValueError: if the file cannot be read, parsed or if the version cannot be determined
    """
    try:
        xml = ElementTree.parse(xml_file)
    except IOError:
        raise ValueError('could not open and or parse the XML file {}'.format(xml_file))

    if is_valid_post_6_2_version(xml):
        return QeXmlVersion.POST_6_2
    elif is_valid_pre_6_2_version(xml):
        return QeXmlVersion.PRE_6_2
    else:
        raise ValueError('unrecognized XML file version')


def get_schema_filepath(xml):
    """
    """
    schema_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schemas')
    schema_filename = get_schema_filename(xml)
    schema_filepath = os.path.join(schema_directory, schema_filename)

    return schema_filepath


def get_default_schema_filepath():
    """
    """
    schema_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schemas')
    schema_filename = DEFAULT_SCHEMA_FILENAME
    schema_filepath = os.path.join(schema_directory, schema_filename)

    return schema_filepath


def get_schema_filename(xml):
    """
    Returns the name of the schema file that corresponds to the given XML object
    """
    
    XML_SCHEMA_LOCATION_KEY = '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'
    # The part in curly brackets is an expanded namespace

    element_root = xml.getroot()

    if element_root is None:
        return None

    element_schema_location = element_root.get(XML_SCHEMA_LOCATION_KEY)
    # e.g. "http://www.quantum-espresso.org/ns/qes/qes-1.0 http://www.quantum-espresso.org/ns/qes/qes-1.0.xsd"

    if element_schema_location is None:
        return None

    schema_location = element_schema_location.split()[1]  # e.g. "http://www.quantum-espresso.org/ns/qes/qes-1.0.xsd"
    schema_filename = schema_location.rpartition('/')[2]  # e.g. "qes-1.0.xsd"
    
    return schema_filename


def is_valid_post_6_2_version(xml):
    """
    Returns whether the given XML object corresponds to an XML output file of Quantum ESPRESSO pw.x post v6.2

    This version of the output XML is parsable with an XSD schema file
    """
    KNOWN_SCHEMA_VERSIONS = [
        'qes-1.0.xsd',
        'qes_181201.xsd',
        'qes_190207.xsd',
        'qes_190304.xsd',  # changes n_opt_steps from positiveInteger to integer
    ]

    schema_filename = get_schema_filename(xml)

    return schema_filename in KNOWN_SCHEMA_VERSIONS


def is_valid_pre_6_2_version(xml):
    """
    Returns whether the given XML object corresponds to an XML output file of Quantum ESPRESSO pw.x pre v6.2
    """
    element_header = xml.find('HEADER')

    if element_header is None:
        return False

    element_format = element_header.find('FORMAT')

    if element_format is None:
        return False

    try:
        name = element_format.attrib['NAME']
    except KeyError:
        return False

    if name != 'QEXML':
        return False

    return True
