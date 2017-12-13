# -*- coding: utf-8 -*-
from collections import Mapping
from aiida.orm.data.parameter import ParameterData

def update_mapping(original, source):
    """
    """
    return_parameter_data = False

    if isinstance(original, ParameterData):
        return_parameter_data = True
        original = original.get_dict()

    if isinstance(source, ParameterData):
        source = source.get_dict()

    for key, value in source.iteritems():
        if (
            key in original 
            and (isinstance(value, Mapping) or isinstance(value, ParameterData))
            and (isinstance(original[key], Mapping) or isinstance(original[key], ParameterData))
        ):
            original[key] = update_mapping(original[key], value)
        else:
            original[key] = value

    if return_parameter_data:
        original = ParameterData(dict=original)

    return original