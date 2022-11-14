# -*- coding: utf-8 -*-
"""Utilities to deal with various mapping data structures."""
from collections.abc import Mapping

from aiida.common import AttributeDict
from aiida.orm import Dict


def get_logging_container():
    """Return an `AttributeDict` that can be used to map logging messages to certain log levels.

    This datastructure is useful to add log messages in a function that does not have access to the right logger. Once
    returned, the caller who does have access to the logger can then easily loop over the contents and pipe the messages
    through the actual logger.

    :return: :py:class:`~aiida.common.extendeddicts.AttributeDict`
    """
    return AttributeDict({
        'debug': [],
        'info': [],
        'warning': [],
        'error': [],
        'critical': [],
    })


def update_mapping(original, source):
    """Update a nested dictionary with another optionally nested dictionary.

    The dictionaries may be plain Mapping objects or `Dict` nodes. If the original dictionary is an instance of `Dict`
    the returned dictionary will also be wrapped in `Dict`.

    :param original: Mapping object or `Dict` instance
    :param source: Mapping object or `Dict` instance
    :return: the original dictionary updated with the source dictionary
    """
    return_node = False

    if isinstance(original, Dict):
        return_node = True
        original = original.get_dict()

    if isinstance(source, Dict):
        source = source.get_dict()

    for key, value in source.items():
        if key in original and isinstance(value, (Dict, Mapping)) and isinstance(original[key], (Dict, Mapping)):
            original[key] = update_mapping(original[key], value)
        else:
            original[key] = value

    if return_node:
        original = Dict(original)

    return original


def prepare_process_inputs(process, inputs):
    """Prepare the inputs for submission for the given process, according to its spec.

    That is to say that when an input is found in the inputs that corresponds to an input port in the spec of the
    process that expects a `Dict`, yet the value in the inputs is a plain dictionary, the value will be wrapped in by
    the `Dict` class to create a valid input.

    :param process: sub class of `Process` for which to prepare the inputs dictionary
    :param inputs: a dictionary of inputs intended for submission of the process
    :return: a dictionary with all bare dictionaries wrapped in `Dict` if dictated by the process spec
    """
    prepared_inputs = wrap_bare_dict_inputs(process.spec().inputs, inputs)
    return AttributeDict(prepared_inputs)


def wrap_bare_dict_inputs(port_namespace, inputs):
    """Wrap bare dictionaries in `inputs` in a `Dict` node if dictated by the corresponding port in given namespace.

    :param port_namespace: a `PortNamespace`
    :param inputs: a dictionary of inputs intended for submission of the process
    :return: a dictionary with all bare dictionaries wrapped in `Dict` if dictated by the port namespace
    """
    from aiida.engine.processes import PortNamespace

    wrapped = {}

    for key, value in inputs.items():

        if key not in port_namespace:
            wrapped[key] = value
            continue

        port = port_namespace[key]
        valid_types = port.valid_type if isinstance(port.valid_type, (list, tuple)) else (port.valid_type,)

        if isinstance(port, PortNamespace):
            wrapped[key] = wrap_bare_dict_inputs(port, value)
        elif Dict in valid_types and isinstance(value, dict):
            wrapped[key] = Dict(value)
        else:
            wrapped[key] = value

    return wrapped
