# -*- coding: utf-8 -*-
"""Utilities for pseudo potentials."""
from __future__ import absolute_import

import six

from aiida.orm.nodes.data.upf import UpfData, get_pseudos_from_structure


def validate_and_prepare_pseudos_inputs(structure, pseudos=None, pseudo_family=None):  # pylint: disable=invalid-name
    """
    Use the explicitly passed pseudos dictionary or use the pseudo_family in combination with
    the structure to obtain that dictionary.

    The pseudos dictionary should now be a dictionary of UPF nodes with the kind as linkname
    As such, if there are multiple kinds with the same element, there will be duplicate UPF nodes
    but multiple links for the same input node are not allowed. Moreover, to couple the UPF nodes
    to the Calculation instance, we have to go through the use_pseudo method, which takes the kind
    name as an additional parameter. When creating a Calculation through a Process instance, one
    cannot call the use methods directly but rather should pass them as keyword arguments. However,
    we can pass the additional parameters by using them as the keys of a dictionary

    :param structure: StructureData node
    :param pseudos: a dictionary where keys are the kind names and value are UpfData nodes
    :param pseudo_family: pseudopotential family name to use, should be Str node
    :raises: ValueError if neither pseudos or pseudo_family is specified or if no UpfData is found for
        every element in the structure
    :returns: a dictionary of UpfData nodes where the key is the kind name
    """
    from aiida.orm import Str

    if pseudos and pseudo_family:
        raise ValueError('you cannot specify both "pseudos" and "pseudo_family"')
    elif pseudos is None and pseudo_family is None:
        raise ValueError('neither an explicit pseudos dictionary nor a pseudo_family was specified')
    elif pseudo_family:
        # This will already raise some exceptions, potentially, like the ones below
        pseudos = get_pseudos_from_structure(structure, pseudo_family.value)
    elif isinstance(pseudos, (six.string_types, Str)):
        raise TypeError('you passed "pseudos" as a string - maybe you wanted to pass it as "pseudo_family" instead?')

    for kind in structure.get_kind_names():
        if kind not in pseudos:
            raise ValueError('no pseudo available for element {}'.format(kind))
        elif not isinstance(pseudos[kind], UpfData):
            raise ValueError('pseudo for element {} is not of type UpfData'.format(kind))

    return pseudos


def get_pseudos_of_calc(calc):
    """
    Return a dictionary of pseudos used by a given (pw.x, cp.x) calculation.

    This returns a dictionary ``pseudos`` that can be set in a builder as ``builder.pseudo = pseudos``.

    :param calc: a pw.x or cp.x calculation.
    :return: a dictionary where the key is the kind name and the value is the UpfData object.
    """
    from aiida.common.links import LinkType

    pseudos = {}
    # I create here a dictionary that associates each kind name to a pseudo
    inputs = calc.get_incoming(link_type=LinkType.INPUT_CALC)
    for linkname in inputs.keys():
        if linkname.startswith(calc._get_linkname_pseudo_prefix()):  # pylint: disable=protected-access
            # Note that this string might be a sequence of kind names
            # concatenated by an underscore, see implementation in the
            # input plugin implementation.
            multiplekindstring = linkname[len(calc._get_linkname_pseudo_prefix()):]  # pylint: disable=protected-access
            pseudos[multiplekindstring] = inputs[linkname]
    return pseudos


def get_pseudos_from_dict(structure, pseudos_uuids):
    """
    Given a dictionary in the format::

        {
            'Al': '045a7a8d-feb1-4aeb-9d32-4e04b13bfc32',
            'C': '08ad7d53-b7cc-45d5-acb8-13530790b751',
        }

    (i.e., associating a chemical element name to a UUID of a UpfData node
    in the database), and a AiiDA structure, return a dictionary
    associating each kind name with its UpfData object.

    :param structure: a StructureData
    :param pseudos_uuids: a dictionary of UUIDs of UpfData for each chemical element, as specified above
    :raise MultipleObjectsError: if more than one UPF for the same element is
       found in the group.
    :raise NotExistent: if no UPF for an element in the group is
       found in the group.
    """
    from aiida.common import NotExistent
    from aiida.orm import load_node

    pseudo_list = {}
    for kind in structure.kinds:
        symbol = kind.symbol
        try:
            uuid = pseudos_uuids[symbol]
        except KeyError:
            raise NotExistent('No UPF for element {} found in the provided pseudos_uuids dictionary'.format(symbol))
        try:
            upf = load_node(uuid)
        except NotExistent:
            raise NotExistent(
                'No node found associated to the UUID {} given for element {} '
                'in the provided pseudos_uuids dictionary'.format(uuid, symbol)
            )
        if not isinstance(upf, UpfData):
            raise ValueError('Node with UUID {} is not a UpfData'.format(uuid))
        if upf.element != symbol:
            raise ValueError(
                'Node<{}> is associated to element {} and not to {} as expected'.format(uuid, upf.element, symbol)
            )

        pseudo_list[kind.name] = upf

    return pseudo_list
