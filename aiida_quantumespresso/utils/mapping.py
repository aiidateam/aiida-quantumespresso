# -*- coding: utf-8 -*-
from collections import Mapping
from aiida.common.extendeddicts import AttributeDict
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


def prepare_process_inputs(inputs):
    """
    Prepare the inputs dictionary for a calculation process. Any remaining bare dictionaries in the inputs
    dictionary will be wrapped in a ParameterData data node except for the '_options' key which should remain
    a standard dictionary. Another exception are dictionaries whose keys are not strings but for example tuples.
    This is the format used by input groups as in for example the explicit pseudo dictionary where the key is
    a tuple of kind to which the UpfData corresponds.
    """
    prepared_inputs = AttributeDict()

    for key, val in inputs.iteritems():
        if key != '_options' and isinstance(val, dict) and all([isinstance(k, (basestring)) for k in val.keys()]):
            prepared_inputs[key] = ParameterData(dict=val)
        else:
            prepared_inputs[key] = val

    return prepared_inputs