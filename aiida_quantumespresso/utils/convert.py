# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function
import numbers
import six
from six.moves import zip


def conv_to_fortran(val, quote_strings=True):
    """
    :param val: the value to be read and converted to a Fortran-friendly string.
    """
    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    import numpy

    if isinstance(val, (bool, numpy.bool_)):
        if val:
            val_str = '.true.'
        else:
            val_str = '.false.'
    elif isinstance(val, numbers.Integral):
        val_str = "{:d}".format(val)
    elif isinstance(val, numbers.Real):
        val_str = ("{:18.10e}".format(val)).replace('e', 'd')
    elif isinstance(val, six.string_types):
        if quote_strings:
            val_str = "'{!s}'".format(val)
        else:
            val_str = "{!s}".format(val)
    else:
        raise ValueError("Invalid value '{}' of type '{}' passed, accepts only bools, ints, floats and strings".format(
            val, type(val)))

    return val_str


def conv_to_fortran_withlists(val, quote_strings=True):
    """
    Same as conv_to_fortran but with extra logic to handle lists
    :param val: the value to be read and converted to a Fortran-friendly string.
    """
    # pylint: disable=too-many-return-statements

    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    if isinstance(val, (list, tuple)):
        val_str = ", ".join(conv_to_fortran(thing, quote_strings=quote_strings) for thing in val)
        return val_str

    if isinstance(val, bool):
        if val:
            return '.true.'

        return '.false.'

    if isinstance(val, six.integer_types):
        return "{:d}".format(val)

    if isinstance(val, float):
        return "{:18.10e}".format(val).replace('e', 'd')

    if isinstance(val, six.string_types):
        if quote_strings:
            return "'{!s}'".format(val)

        return "{!s}".format(val)

    raise ValueError("Invalid value passed, accepts only bools, ints, floats and strings")


def convert_input_to_namelist_entry(key, val, mapping=None):
    """
    Given a key and a value, from an input parameters dictionary for a namelist calculation, map it to
    the appropriate string format for the namelist input file. For single values it will return a single
    string, but for values that are a dictionary, list or tuple, the returned string may be multiline.

    :param key: the namelist keyword name
    :param val: the namelist keyword value
        The value can be either a single value, list/tuple, a double nested list or a dictionary.
        Depending on the type of the value the resulting string differs vastly

        * single list:
            A list of keywords will be generated, where the index of the value in the list will be
            used as the index in the keyword and the value itself will be converted using conv_to_fortran.
            For example::

                'efield': [4, 5, 6]

            will result in::

                efield(1) = 4
                efield(1) = 5
                efield(1) = 6

        * double nested list:
            This format can be used for keywords that require one or more indices that do not necessarily
            follow a sequential number, but take specific values that need to be defined by the user.
            For example::

                'starting_ns_eigenvalue': [
                    [1, 1, 3, 3.5],
                    [2, 1, 1, 2.8]
                ]

            will be formatted as::

                starting_ns_eigenvalue(1,1,3) = 3.5
                starting_ns_eigenvalue(2,1,1) = 2.8

            Note that if the mapping argument is provided in the input, any value in sub lists that matches
            a key in the mapping dictionary (that is to say it is a string that matches one of the kinds), it
            will be replaced with the index of the corresponding atomic species. For example::

                hubbard_j: [
                    [2, 'Ni', 3.5],
                    [2, 'Fe', 7.4],
                ]

            would be formatted as::

                hubbard_j(2, 1) = 3.5
                hubbard_j(2, 3) = 7.4

            Assuming the mapping dictionary contained the kinds 'Ni' and 'Fe', with the indices 1 and 3, respectively

        * dictionary:
            The keys of this dictionary should correspond to keys in the mapping input argument and will be replaced
            with the corresponding value. This can be used for keywords that take a single index that needs to conform
            to the index of the atomic species to which the keyword value should apply. For example::

                hubbard_u: {
                    'Co': 3.5,
                    'O': 7.4,
                }

            will be formatted as::

                hubbard_u(1) = 3.5
                hubbard_u(3) = 7.4

            assuming that the kinds 'Co' and 'O' would have atomic species indices 1 and 3, respectively.
            This mapping from kind name to atomic species index should be defined by the `mapping` argument.

    :param mapping: optional parameter, that must be provided if val is a dictionary or a double nested list
        where the sub lists contain string values. The keys of the mapping dictionary should be the atomic species
        names that will be encountered in the value, and the corresponding value should be the index of that
        atomic species. Example::

            mapping = {
                'Fe': 1,
                'O': 2,
            }

        This will map every occurrence of 'Fe' and 'O' in the values to the corresponding integer.
    """
    # I don't try to do iterator=iter(val) and catch TypeError because it would also match strings
    # I check first the dictionary, because it would also match hasattr(__iter__)
    if isinstance(val, dict):

        if mapping is None:
            raise ValueError("If 'val' is a dictionary, you must provide also the 'mapping' parameter")

        # At difference with the case of a list, at the beginning list_of_strings
        # is a list of 2-tuples where the first element is the idx, and the
        # second is the actual line. This is used at the end to resort everything.
        list_of_strings = []

        for elemk, itemval in six.iteritems(val):
            try:
                idx = mapping[elemk]
            except KeyError:
                raise ValueError("Unable to find the key '{}' in the mapping dictionary".format(elemk))

            list_of_strings.append((idx, u'  {0}({2}) = {1}\n'.format(key, conv_to_fortran(itemval), idx)))

        # I first have to resort, then to remove the index from the first column, finally to join the strings
        list_of_strings = list(zip(*sorted(list_of_strings)))[1]
        return ''.join(list_of_strings)

    # A list/tuple of values
    elif isinstance(val, (list, tuple)):

        list_of_strings = []

        for idx, itemval in enumerate(val):

            if isinstance(itemval, (list, tuple)):

                values = []

                for value in itemval[:-1]:

                    if not isinstance(value, (int, six.string_types)):
                        raise ValueError('values of double nested lists should be either integers or strings')

                    if isinstance(value, six.string_types):
                        if mapping is None:
                            raise ValueError('cannot map the string value because no mapping was defined')

                        if value not in mapping:
                            raise ValueError('the nested list contained string {} but this is not a key in the mapping'
                                .format(value))
                        else:
                            values.append(str(mapping[value]))
                    else:
                        values.append(str(value))

                idx_string = ','.join(values)
                itemval = itemval.pop()
            else:
                idx_string = '{}'.format(idx + 1)

            list_of_strings.append(u'  {0}({2}) = {1}\n'.format(key, conv_to_fortran(itemval), idx_string))

        return u''.join(list_of_strings)

    # Single value
    else:
        return u'  {0} = {1}\n'.format(key, conv_to_fortran(val))
