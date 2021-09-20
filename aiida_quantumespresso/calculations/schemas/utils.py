# -*- coding: utf-8 -*-
"""Utilities for providing JSON schemas for the validation of calculation inputs."""

import itertools
import re
import random

XML_TYPE_MAPPING = {
    'json': {
        'CHARACTER': 'string',
        'LOGICAL': 'boolean',
        'INTEGER': 'integer',
        'REAL': 'number',
    },
}


def convert_to_json_schema(input_xml, schema_version=None):
    """Convert a Quantum ESPRESSO input XML file into a JSON schema."""

    def _convert_type(qe_type, to_type='json'):
        """Convert a QE type in the input XML file to a JSON schema type."""
        return XML_TYPE_MAPPING[to_type][qe_type.upper()]

    def _get_allowed_values(variable, var_type):
        """Get the allowed values of a variable in the xml schema."""

        allowed_values = [opt.getAttribute('val').split(', ') for opt in variable.getElementsByTagName('opt')]
        if var_type == 'integer':
            allowed_values = [int(re.sub('[^0-9]', '', value)) for value in itertools.chain(*allowed_values)]
        else:
            allowed_values = [value.strip('\' ') for value in itertools.chain(*allowed_values)]

        return allowed_values

    schema_version = schema_version or 'http://json-schema.org/draft/2019-09/schema#'

    json_schema = {'$schema': schema_version, 'type': 'object', 'additionalProperties': False, 'properties': {}}

    for namelist in input_xml.getElementsByTagName('namelist'):

        namelist_dict = {'type': 'object', 'additionalProperties': False, 'properties': {}}
        for tagname in ('var', 'dimension', 'multidimension'):

            for var in namelist.getElementsByTagName(tagname):

                var_dict = {}

                if tagname == 'dimension':
                    var_dict['type'] = 'object'  # These are provided as a dict in the plugin
                elif tagname == 'multidimension':
                    var_dict['type'] = 'array'  # These are provided as a list in the plugin
                    var_dict['items'] = {'type': 'array', 'items': {'anyOf': [{'type': 'number'}, {'type': 'string'}]}}
                else:
                    try:
                        var_dict['type'] = _convert_type(var.getAttribute('type'))
                    except KeyError:
                        continue  # Skip variables without a type -> These are stored in vargoups

                allowed_values = _get_allowed_values(var, var_dict['type'])

                if len(allowed_values) > 0:
                    var_dict['enum'] = allowed_values

                namelist_dict['properties'][var.getAttribute('name')] = var_dict

        for vargroup in namelist.getElementsByTagName('vargroup'):

            for var in vargroup.getElementsByTagName('var'):

                var_dict = {'type': _convert_type(vargroup.getAttribute('type'))}
                namelist_dict['properties'][var.getAttribute('name')] = var_dict

        json_schema['properties'][namelist.getAttribute('name')] = namelist_dict

    return json_schema


def generate_params_from_schema(json_schema):
    """Generate a QE ``parameter`` dictionary from a JSON schema."""

    parameters = {}

    for namelist in ('CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'):

        parameters[namelist] = {}

        for tag, schema in json_schema['properties'][namelist]['properties'].items():

            if 'enum' in schema:
                parameters[namelist][tag] = random.choice(schema['enum'])
                continue
            if schema['type'] == 'object':
                parameters[namelist][tag] = {'a': 1}
            if schema['type'] == 'array':
                parameters[namelist][tag] = [[0, 1], [2, 3]]
            if schema['type'] == 'string':
                parameters[namelist][tag] = 'Geef mij maar zo`n lekker stringetje'
            if schema['type'] == 'integer':
                parameters[namelist][tag] = 7
            if schema['type'] == 'boolean':
                parameters[namelist][tag] = False
            if schema['type'] == 'number':
                parameters[namelist][tag] = 3.14

    return parameters
