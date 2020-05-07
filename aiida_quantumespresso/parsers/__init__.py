# -*- coding: utf-8 -*-
from aiida.common import OutputParsingError


class QEOutputParsingError(OutputParsingError):
    """Exception raised when there is a parsing error in the QE parser."""
    pass


def get_parser_info(parser_info_template=None):
    """Return a template dictionary with details about the parser such as the version.

    :param parser_info_template: template string with single placeholder to be replaced by current version number
    :returns: dictionary with parser name, version and empty list for warnings
    """
    import aiida_quantumespresso

    parser_version = aiida_quantumespresso.__version__
    parser_info = {}
    parser_info['parser_warnings'] = []
    parser_info['parser_version'] = parser_version

    if parser_info_template is None:
        parser_info['parser_info'] = 'aiida-quantumespresso parser v{}'.format(parser_version)
    else:
        parser_info['parser_info'] = parser_info_template.format(parser_version)

    return parser_info
