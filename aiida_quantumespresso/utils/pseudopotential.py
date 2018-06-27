# -*- coding: utf-8 -*-
from aiida.orm.data.upf import UpfData, get_pseudos_from_structure


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
    inputs = calc.get_inputs_dict(link_type=LinkType.INPUT)
    for linkname in inputs.keys():
        if linkname.startswith(calc._get_linkname_pseudo_prefix()):
            # Note that this string might be a sequence of kind names
            # concatenated by an underscore, see implementation in the
            # input plugin implementation.
            multiplekindstring = linkname[len(calc._get_linkname_pseudo_prefix()):]
            pseudos[multiplekindstring] = inputs[linkname]
    return pseudos


def validate_and_prepare_pseudos_inputs(structure, pseudos=None, pseudo_family=None):
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
    :param pseudo_family: string name of the pseudopotential family to use
    :raises: ValueError if neither pseudos or pseudo_family is specified or if no UpfData is found for
        every element in the structure
    :returns: a dictionary of UpfData nodes where the key is a tuple with the kind name
    """
    from aiida.orm.data.base import Str
    result_pseudos = {}
    unique_pseudos = {}

    if pseudos and pseudo_family:
        raise ValueError('You cannot specify both "pseudos" and "pseudo_family"')
    elif pseudos is None and pseudo_family is None:
        raise ValueError('Neither an explicit pseudos dictionary nor a pseudo_family was specified')
    elif pseudo_family:
        # This will already raise some exceptions, potentially, like the ones below
        pseudos = get_pseudos_from_structure(structure, str(pseudo_family))

    if isinstance(pseudos, (str, unicode, Str)):
        raise TypeError('You passed "pseudos" as a string - maybe you wanted to pass it as "pseudo_family" instead?')

    for kind in structure.get_kind_names():
        if kind not in pseudos:
            raise ValueError('no pseudo available for element {}'.format(kind))
        elif not isinstance(pseudos[kind], UpfData):
            raise ValueError('pseudo for element {} is not of type UpfData'.format(kind))

    for kind, pseudo in pseudos.iteritems():
        unique_pseudos.setdefault(pseudo, []).append(kind)

    for pseudo, kinds in unique_pseudos.iteritems():
        result_pseudos[tuple(kinds)] = pseudo

    return result_pseudos
