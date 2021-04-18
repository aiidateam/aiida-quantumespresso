# -*- coding: utf-8 -*-
"""Utilities for providing JSON schemas for the validation of calculation inputs."""

import itertools

def convert_to_json_schema(input_xml, schema_version=None):
    """Convert a Quantum ESPRESSO input XML file into a JSON schema."""
    
    def _convert_type(qe_type):
        """Convert a QE type in the input XML file to a JSON schema type."""
        type_mapping = {
            'CHARACTER': 'string',
            'LOGICAL': 'boolean',
            'INTEGER': 'integer',
            'REAL': 'number',
        }
        return type_mapping[qe_type.upper()]

    schema_version = schema_version or 'http://json-schema.org/draft/2019-09/schema#'
    
    json_schema = {
        '$schema': schema_version,
        'type': 'object',
        'additionalProperties': False,
        'properties': {}
    }
    
    for namelist in input_xml.getElementsByTagName('namelist'):

        namelist_dict = {
            'type': 'object',
            'additionalProperties': False,
            'properties': {}
        }
        for tagname in ('var', 'dimension', 'multidimension'):
            
            for var in namelist.getElementsByTagName(tagname):
                
                if tagname == 'dimension':
                    var_type = 'object'  # These are provided as a dict in the plugin
                elif tagname == 'multidimension':
                    var_type = 'array'  # These are provided as a list in the plugin
                else:
                    try:
                        var_type = _convert_type(var.getAttribute('type'))
                    except KeyError:
                        continue  # Skip variables without a type -> These are stored in vargoups
                allowed_values = [
                    opt.getAttribute('val').split(', ') for opt in var.getElementsByTagName('opt')
                ]
                allowed_values = [value.strip('\' ') for value in itertools.chain(*allowed_values)]
                
                var_dict = {
                    'type': var_type,
                }
                if len(allowed_values) > 0:
                    var_dict['enum'] = allowed_values
                
                namelist_dict['properties'][var.getAttribute('name')] = var_dict

        for vargroup in namelist.getElementsByTagName('vargroup'):
            
            for var in vargroup.getElementsByTagName('var'):
            
                var_dict = {
                    'type': _convert_type(vargroup.getAttribute('type'))
                }
                namelist_dict['properties'][var.getAttribute('name')] = var_dict
            
        json_schema['properties'][namelist.getAttribute('name')] = namelist_dict

    return json_schema
