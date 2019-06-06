# -*- coding: utf-8 -*-
from __future__ import absolute_import
from collections import Mapping
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.nodes.data.dict import Dict
import six


def update_mapping(original, source):
    """
    Update a nested dictionary with another optionally nested dictionary. The dictionaries may be plain
    Mapping objects or Dict nodes. If the original dictionary is an instance of Dict
    the returned dictionary will also be wrapped in Dict.

    :param original: Mapping object or Dict instance
    :param source: Mapping object or Dict instance
    :return: the original dictionary updated with the source dictionary
    """
    return_parameter_data = False

    if isinstance(original, Dict):
        return_parameter_data = True
        original = original.get_dict()

    if isinstance(source, Dict):
        source = source.get_dict()

    for key, value in six.iteritems(source):
        if (
            key in original and
            (isinstance(value, Mapping) or isinstance(value, Dict)) and
            (isinstance(original[key], Mapping) or isinstance(original[key], Dict))
        ):
            original[key] = update_mapping(original[key], value)
        else:
            original[key] = value

    if return_parameter_data:
        original = Dict(dict=original)

    return original


def prepare_process_inputs(process, inputs):
    """
    Prepare the inputs for submission for the given process, according to its spec. That is to say that
    when an input is found in the inputs that corresponds to an input port in the spec of the process that
    expects a Dict, yet the value in the inputs is a plain dictionary, the value will be wrapped
    in by the Dict class to create a valid input.

    :param process: sub class of Process for which to prepare the inputs dictionary
    :param inputs: a dictionary of inputs intended for submission of the process
    :return: a dictionary with all bare dictionaries wrapped in Dict if dictated by process spec
    """
    prepared_inputs = AttributeDict()

    try:
        process_spec = process.spec()
    except AttributeError:
        raise ValueError('process {} does not have a spec')

    for key, value in six.iteritems(inputs):

        if key not in process_spec.inputs:
            continue

        if process_spec.inputs[key].valid_type == Dict and isinstance(value, dict):
            prepared_inputs[key] = Dict(dict=value)
        else:
            prepared_inputs[key] = value

    return prepared_inputs
