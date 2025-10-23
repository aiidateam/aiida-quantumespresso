"""Wrapper module for routing XML parsing to legacy or modern parsers based on QE version.

.. deprecated:: 4.6
    This module has been deprecated and will be removed in aiida-quantumespresso v5.0.
    Use parse_xml() from parse_xml.parse directly. Legacy XML format support (QE < v6.2) will be dropped.
"""

import warnings
from xml.etree import ElementTree

from aiida.common.warnings import AiidaDeprecationWarning

from aiida_quantumespresso.parsers.parse_xml.exceptions import XMLParseError
from aiida_quantumespresso.parsers.parse_xml.parse import parse_xml_post_6_2
from aiida_quantumespresso.parsers.parse_xml.versions import QeXmlVersion, get_xml_file_version

from .legacy import parse_pw_xml_pre_6_2

warnings.warn(
    'The parse_xml.pw.parse module has been deprecated and will be removed in aiida-quantumespresso v5.0.\n'
    'Use parse_xml() from parse_xml.parse directly. Legacy XML format support will be dropped.',
    AiidaDeprecationWarning,
)


def parse_xml(xml_file, dir_with_bands=None):
    try:
        xml_parsed = ElementTree.parse(xml_file)
    except ElementTree.ParseError as exception:
        raise XMLParseError('error while parsing XML file') from exception

    xml_file_version = get_xml_file_version(xml_parsed)

    if xml_file_version == QeXmlVersion.POST_6_2:
        parsed_data, logs = parse_xml_post_6_2(xml_parsed)
    elif xml_file_version == QeXmlVersion.PRE_6_2:
        xml_file.seek(0)
        parsed_data, logs = parse_pw_xml_pre_6_2(xml_file, dir_with_bands)

    return parsed_data, logs
