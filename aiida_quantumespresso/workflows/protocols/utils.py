# -*- coding: utf-8 -*-
"""Utilities to manipulate the workflow input protocols."""
import functools
import os
import pathlib
import yaml


def recursive_merge(left, right):
    """Recursively merge two dictionaries into a single dictionary.

    If any key is present in both ``left`` and ``right`` dictionaries, the value from the ``right`` dictionary is
    assigned to the key.

    :param left: first dictionary
    :param right: second dictionary
    :return: the recursively merged dictionary
    """
    import collections

    for key, value in left.items():
        if key in right:
            if isinstance(value, collections.abc.Mapping) and isinstance(right[key], collections.abc.Mapping):
                right[key] = recursive_merge(value, right[key])

    merged = left.copy()
    merged.update(right)

    return merged


def get_protocol_inputs(cls, protocol=None, overrides=None):
    """Docs."""
    from aiida.plugins.entry_point import get_entry_point_from_class

    _, entry_point = get_entry_point_from_class(cls.__module__, cls.__name__)
    entry_point_name = entry_point.name
    parts = entry_point_name.split('.')
    parts.pop(0)
    filename = f'{parts.pop()}.yaml'
    basepath = functools.reduce(os.path.join, parts)

    with (pathlib.Path(__file__).resolve().parent / basepath / filename).open() as handle:
        data = yaml.safe_load(handle)

    protocol = protocol or data['default_protocol']
    inputs = recursive_merge(data['default_inputs'], data['protocols'][protocol])
    inputs.pop('description')

    if overrides:
        return recursive_merge(inputs, overrides)

    return inputs
