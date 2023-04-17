# -*- coding: utf-8 -*-
import enum
import os

from aiida_quantumespresso.parsers.parse_xml.exceptions import XMLUnsupportedFormatError

DEFAULT_SCHEMA_FILENAME = 'qes-1.0.xsd'
DIRNAME_SCHEMAS = 'schemas'
DIRPATH_SCHEMAS = os.path.join(os.path.dirname(os.path.abspath(__file__)), DIRNAME_SCHEMAS)


class QeXmlVersion(enum.Enum):
    """An enum with the versions of XML output file known to exist for Quantum ESPRESSO."""

    PRE_6_2 = 0
    POST_6_2 = 1


def get_xml_file_version(xml):
    """Return the version of the Quantum ESPRESSO pw.x and cp.x XML output file.

    :param xml: the pre-parsed XML object
    :raises XMLUnsupportedFormatError: if the file cannot be read, parsed or if the version cannot be determined
    """
    if is_valid_post_6_2_version(xml):
        return QeXmlVersion.POST_6_2
    elif is_valid_pre_6_2_version(xml):
        return QeXmlVersion.PRE_6_2
    else:
        raise XMLUnsupportedFormatError(f'unrecognized XML file version: cannot find schema {get_schema_filename(xml)} in {os.path.join(os.path.dirname(os.path.abspath(__file__)), "schemas")}. You can look for it in https://github.com/QEF/qeschemas')


def get_schema_filepath(xml):
    """Return the absolute filepath to the XML schema file that can be used to parse the given XML.

    :param xml: the pre-parsed XML object
    :return: the XSD absolute filepath
    """
    schema_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schemas')
    schema_filename = get_schema_filename(xml)
    schema_filepath = os.path.join(schema_directory, schema_filename)

    return schema_filepath


def get_default_schema_filepath():
    """Return the absolute filepath to the default XML schema file.

    :param xml: the pre-parsed XML object
    :return: the XSD absolute filepath
    """
    schema_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schemas')
    schema_filename = DEFAULT_SCHEMA_FILENAME
    schema_filepath = os.path.join(schema_directory, schema_filename)

    return schema_filepath


def get_schema_filename(xml):
    """Return the filename of the schema file that corresponds to the given XML object.

    :param xml: the pre-parsed XML object
    :return: the XSD filename
    """
    xml_schema_location_key = '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'
    # The part in curly brackets is an expanded namespace

    element_root = xml.getroot()

    if element_root is None:
        return None

    # Patch for QE v7.0: The scheme file name was not updated in the `xsi.schemaLocation` element
    try:
        if element_root.find('general_info').find('creator').get('VERSION') == '7.0':
            return 'qes_211101.xsd'
    except AttributeError:
        pass

    element_schema_location = element_root.get(xml_schema_location_key)
    # e.g. "http://www.quantum-espresso.org/ns/qes/qes-1.0 http://www.quantum-espresso.org/ns/qes/qes-1.0.xsd"

    if element_schema_location is None:
        return None

    schema_location = element_schema_location.split()[1]  # e.g. "http://www.quantum-espresso.org/ns/qes/qes-1.0.xsd"
    schema_filename = schema_location.rpartition('/')[2]  # e.g. "qes-1.0.xsd"

    return schema_filename


def get_available_xml_schemas():
    """Return the filenames of the available XML schemas.

    These are essentially all the files with the .xsd extenstion in the `DIRPATH_SCHEMAS` folder

    :return: a list of XML schema filenames
    """
    return [file for file in os.listdir(DIRPATH_SCHEMAS) if file.endswith('.xsd')]


def is_valid_post_6_2_version(xml):
    """Return whether the given XML object corresponds to an XML output file of Quantum ESPRESSO pw.x post v6.2.

    These versions of the output XML are parsable with an XSD schema file

    :param xml: a parsed XML output file
    :return: boolean, True when the XML is parsable with an XSD schema file, False otherwise
    """
    return get_schema_filename(xml) in get_available_xml_schemas()


def is_valid_pre_6_2_version(xml):
    """Returns whether the given XML object corresponds to an XML output file of Quantum ESPRESSO pw.x pre v6.2.

    :param xml: a parsed XML output file
    :return: boolean, True when the XML was produced by Quantum ESPRESSO with the old XML format
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
