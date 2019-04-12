# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.common import EntryPointError

try:
    from reentry import manager as epm
except ImportError:
    import pkg_resources as epm


def get_plugins(category):
    """
    Get a list of plugins for the given category
    """
    return [ep.name for ep in epm.iter_entry_points(group=category)]


def get_plugin(category, name):
    """
    Return an instance of the class registered under the given name and
    for the specified plugin category.

    :param category: the plugin category to load the plugin from, e.g. 'transports'.
    :param name: the name of the plugin
    """
    eps = [ep for ep in epm.iter_entry_points(group=category) if ep.name == name]

    if not eps:
        raise EntryPointError(
            "No plugin named '{}' found for '{}'".format(name, category))

    if len(eps) > 1:
        raise EntryPointError(
            "Multiple plugins found for '{}' in '{}'".format(name, category))

    entrypoint = eps[0]

    try:
        plugin = entrypoint.load()
    except ImportError:
        import traceback
        raise EntryPointError("Loading the plugin '{}' failed:\n{}"
            .format(name, traceback.format_exc()))

    return plugin
