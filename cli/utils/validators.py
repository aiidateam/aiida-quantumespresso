# -*- coding: utf-8 -*-
from __future__ import absolute_import 
import click

def validate_structure(ctx, param, value):
    """
    Command line option validator for an AiiDA structure data pk. It expects
    an integer for the value and will try to load the corresponding node. it
    will also check if successful if the node is a StructureData instance.

    :param value: a StructureData node pk
    :returns: a StructureData instance
    """
    from aiida.common.exceptions import NotExistent
    from aiida.orm import load_node
    from aiida.orm.data.structure import StructureData

    try:
        structure = load_node(int(value))
    except NotExistent as exception:
        raise click.BadParameter('failed to load the node<{}>\n{}'.format(value, exception))

    if not isinstance(structure, StructureData):
        raise click.BadParameter('node<{}> is not of type StructureData'.format(value))

    return structure

def validate_code(ctx, param, value):
    """
    Command line option validator for an AiiDA code. It expects a string for the value
    that corresponds to the label of an AiiDA code that has been setup.

    :param value: a Code label
    :returns: a Code instance
    """
    from aiida.common.exceptions import NotExistent
    from aiida.orm import Code

    try:
        code = Code.get_from_string(value)
    except NotExistent as exception:
        raise click.BadParameter("failed to load the code with the label '{}'\n{}".format(value, exception))

    return code

def validate_pseudo_family(ctx, param, value):
    """
    Command line option validator for a pseudo potential family. The value should be a
    string that corresponds to a registered UpfData family, which is an AiiDA Group.

    :param value: a UpfData pseudo potential family label
    :returns: the pseudo potential family label
    """
    from aiida.common.exceptions import NotExistent
    from aiida.orm.data.upf import UpfData

    try:
        pseudo_family = UpfData.get_upf_group(value)
    except NotExistent as exception:
        raise click.BadParameter("failed to load the pseudo family the label '{}'\n{}".format(value, exception))

    return value

def validate_kpoint_mesh(ctx, param, value):
    """
    Command line option validator for a kpoints mesh tuple. The value should be a tuple
    of three positive integers out of which a KpointsData object will be created with
    a mesh equal to the tuple.

    :param value: a tuple of three positive integers
    :returns: a KpointsData instance
    """
    from aiida.orm.data.array.kpoints import KpointsData

    if any([type(i) != int for i in value]) or any([int(i) <= 0 for i in value]):
        raise click.BadParameter('all values of the tuple should be positive greater than zero integers')

    try:
        kpoints = KpointsData()
        kpoints.set_kpoints_mesh(value)
    except ValueError as exception:
        raise click.BadParameter("failed to create a KpointsData mesh out of {}\n{}".format(value, exception))

    return kpoints