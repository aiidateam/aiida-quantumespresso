# -*- coding: utf-8 -*-
"""Utilities to automatically format, convert and validate data structures from python to fortran."""
import copy
import difflib
import os
import xml.dom.minidom

from aiida.common import InputValidationError, InternalError
from packaging.version import Version


class QEInputValidationError(InputValidationError):
    """Raise when the parser encounters an error while creating the input file of Quantum ESPRESSO."""


def _check_and_convert(keyword, val, expected_type):
    """Check and convert the value for the given keyword against an expected type.

    :param keyword: the keyword
    :param val: value for the given keyword
    :param expected_type: type the value should have: [`INTEGER`, `REAL`, `CHARACTER`, `LOGICAL`]
    """
    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    if expected_type.upper() == 'LOGICAL':
        if isinstance(val, bool):
            outval = val
        else:
            raise TypeError(f'Expected a boolean for keyword {keyword}, found {type(val)} instead')
    elif expected_type.upper() == 'REAL':
        if isinstance(val, int):
            outval = float(val)
        elif isinstance(val, float):
            outval = val
        else:
            raise TypeError(f'Expected a float for keyword {keyword}, found {type(val)} instead')
    elif expected_type.upper() == 'INTEGER':
        if isinstance(val, int):
            outval = val
        else:
            raise TypeError(f'Expected an integer for keyword {keyword}, found {type(val)} instead')
    elif expected_type.upper() == 'CHARACTER':
        if isinstance(val, str):
            outval = val
        else:
            raise TypeError(f'Expected a string for keyword {keyword}, found {type(val)} instead')
    else:
        raise InternalError(f'Unexpected type check for keyword {keyword}: {expected_type.upper()})')

    return outval


def pw_input_helper(input_params, structure, stop_at_first_error=False, flat_mode=False, version='6.2'):
    """Validate if the input dictionary for Quantum ESPRESSO is valid.

    Return the dictionary (possibly with small variations: e.g. convert integer to float where necessary, recreate the
    proper structure if flat_mode is True, ...) to use as input parameters (use_parameters) for the pw.x calculation.

    :param input_params: If flat_mode is True, pass a dictionary
        with 'key' = value; use the correct type
        (int, bool, ...) for value. If an array is required:

           * if its length is fixed: pass a list of the required length

           * if its length is 'ntyp': pass a dictionary, associating each
             specie to its value.

           * (other lengths are not supported)

       Example::

             {
             'calculation': 'vc-relax',
             'ecutwfc': 30.,
             'hubbard_u': {'O': 1},
             }

       If instead flat_mode is False, pass a dictionary in the format
       expected by AiiDA (keys are namelists, values are in the format
       specified above, i.e. key/value pairs for all keywords in the
       given namelist).
       Example::

             {
                 'CONTROL': {
                     'calculation': 'vc-relax'
                     },
                 'SYSTEM': {
                     'hubbard_u': {'O': 1.0},
                     'ecutwfc': 30.,
                     },
             },


    :param structure: the StructureData object used as input for QE pw.x
    :param stop_at_first_error: if True, stops at the first error.
        Otherwise, when, possible, continue and give a global error for all
        the issues encountered.
    :param flat_mode: if True, instead of passing the dictionary of namelists,
        and inside the keywords, pass directly the keywords - this function
        will return the correct dictionary to pass to the PwCalculation,
        with the keywords arranged in the correct namelist.
    :param version: string with version number, used to find the correct XML
        file descriptor. If not specified, uses the most recent version
        available in the validator. It reads the definitions from the XML files
        in the same folder as this python module. If the version is not
        recognised, the Exception message will also suggest a close-by version.

    :raise QEInputValidationError:
        if the input is not considered valid.
    """
    # pylint: disable=too-many-branches,too-many-statements
    errors_list = []

    # =========== LIST OF KNOWN NAMELISTS, CARDS, VARIABLES, ... ===============
    compulsory_namelists = ['CONTROL', 'SYSTEM', 'ELECTRONS']

    valid_calculations_and_opt_namelists = {  # pylint: disable=invalid-name
        'scf': [],
        'nscf': [],
        'bands': [],
        'relax': ['IONS'],
        'md': ['IONS'],
        'vc-relax': ['IONS', 'CELL'],
        'vc-md': ['IONS', 'CELL'],
    }

    if not isinstance(input_params, dict):
        raise QEInputValidationError('input_params must be a dictionary')
    # So that if I modify input_params, nothing happens outside
    if flat_mode:
        input_params_internal = copy.deepcopy(input_params)
    else:
        input_params_internal = {}
        input_original_namelists = {}
        all_input_namelists = set()
        for namelist, content in input_params.items():
            if not isinstance(content, dict):
                raise QEInputValidationError(
                    f"The content associated to the namelist '{namelist}' must be a dictionary"
                )
            all_input_namelists.add(namelist)
            for key, value in content.items():
                input_params_internal[key] = copy.deepcopy(value)
                if key in input_original_namelists:
                    err_str = f"Keyword '{key}' was specified both in {input_original_namelists[key]} and {namelist}."
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                input_original_namelists[key] = namelist

    # List of the keywords that must not appear in the input
    # (e.g. because they are automatically filled in by the plugin)
    blocked_kws = [
        i.lower() for i in [
            'pseudo_dir',
            'outdir',
            'celldm',
            'nat',
            'ntyp',
            'prefix',
            'a',
            'b',
            'c',
            'cosab',
            'cosac',
            'cosbc',
        ]
    ]

    # List of the keywords that must ALWAYS appear in the input
    compulsory_kws = {i.lower() for i in [
        'calculation',
        'ecutwfc',
    ]}

    # ===================== PARSING OF THE XML DEFINITION FILE ===============
    module_dir = os.path.dirname(__file__)
    if module_dir == '':
        module_dir = os.curdir
    xml_path = os.path.join(module_dir, f'INPUT_PW-{version}.xml')
    try:
        with open(xml_path, 'r', encoding='utf-8') as handle:
            dom = xml.dom.minidom.parse(handle)
    except IOError as exception:
        prefix = 'INPUT_PW-'
        suffix = '.xml'
        versions = [
            fname[len(prefix):-len(suffix)]
            for fname in os.listdir(module_dir)
            if fname.startswith(prefix) and fname.endswith(suffix)
        ]
        versions = sorted(versions, key=Version)
        strictversions = versions + [version]
        strictversions = sorted(strictversions, key=Version)
        pos = strictversions.index(version)
        if pos == 0:
            add_str = ' (the version you specified is too old)'
        else:
            add_str = f' (the older, closest version you can use is {strictversions[pos - 1]})'
        raise QEInputValidationError(
            f"Unknown Quantum Espresso version: {version}. Available versions: {', '.join(versions)};{add_str}"
        ) from exception

    # ========== List of known PW variables (from XML file) ===============
    known_kws = dom.getElementsByTagName('var')
    valid_kws = {}
    for keyword in known_kws:
        if keyword in valid_kws:
            raise InternalError(
                f"Something strange, I found more than one keyword '{keyword}' in the XML description..."
            )

        valid_kws[keyword.getAttribute('name').lower()] = {}
        parent = keyword
        try:
            while True:
                parent = parent.parentNode
                if parent.tagName == 'namelist':
                    valid_kws[keyword.getAttribute('name').lower()]['namelist'] = parent.getAttribute('name').upper()
                    break
        except AttributeError:
            # There are also variables in cards instead of namelists:
            # I ignore them
            pass
            # raise QEInputValidationError("Unable to find namelist for "
            #     "keyword %s." % kw.getAttribute('name'))
        expected_type = keyword.getAttribute('type')
        # Fix for groups of variables
        if expected_type == '':
            if keyword.parentNode.tagName == 'vargroup':
                expected_type = keyword.parentNode.getAttribute('type')
        valid_kws[keyword.getAttribute('name').lower()]['expected_type'] = expected_type.upper()

    # ====== List of known PW 'dimensions' (arrays) (from XML file) ===========
    known_dims = dom.getElementsByTagName('dimension')
    valid_dims = {}
    for dim in known_dims:
        if dim in valid_dims:
            raise InternalError(f"Something strange, I found more than one keyword '{dim}' in the XML description...")

        valid_dims[dim.getAttribute('name').lower()] = {}
        parent = dim
        try:
            while True:
                parent = parent.parentNode
                if parent.tagName == 'namelist':
                    valid_dims[dim.getAttribute('name').lower()]['namelist'] = parent.getAttribute('name').upper()
                    break
        except AttributeError:
            # There are also variables in cards instead of namelists:
            # I ignore them
            pass
            # raise QEInputValidationError("Unable to find namelist "
            #     "for keyword %s." % dim.getAttribute('name'))
        expected_type = dim.getAttribute('type')
        # Fix for groups of variables
        if expected_type == '':
            if dim.parentNode.tagName == 'vargroup':
                expected_type = dim.parentNode.getAttribute('type')
        valid_dims[dim.getAttribute('name').lower()]['expected_type'] = expected_type.upper()
        # I assume start_val is always 1
        start_val = dim.getAttribute('start')
        if start_val != '1':
            raise InternalError(
                f"Wrong start value '{start_val}' in input array (dimension) {dim.getAttribute('name')}"
            )
        # I save the string as it is; somewhere else I will check for its value
        valid_dims[dim.getAttribute('name').lower()]['end_val'] = dim.getAttribute('end')

    # ====== List of known PW 'multidimensions' (arrays) (from XML file) ===========
    known_multidims = dom.getElementsByTagName('multidimension')
    valid_multidims = {}
    for dim in known_multidims:
        if dim in valid_multidims:
            raise InternalError(
                f"Something strange, I found more than one multidimensional keyword '{dim}' in the XML description..."
            )

        valid_multidims[dim.getAttribute('name').lower()] = {}
        parent = dim
        try:
            while True:
                parent = parent.parentNode
                if parent.tagName == 'namelist':
                    valid_multidims[dim.getAttribute('name').lower()]['namelist'] = parent.getAttribute('name').upper()
                    break
        except AttributeError:
            # There are also variables in cards instead of namelists: ignore them
            pass

        expected_type = dim.getAttribute('type').upper()
        start_values = dim.getAttribute('start').split(',')
        end_values = dim.getAttribute('end').split(',')
        indexes = dim.getAttribute('indexes').split(',')

        valid_multidims[dim.getAttribute('name').lower()]['expected_type'] = expected_type
        valid_multidims[dim.getAttribute('name').lower()]['start'] = start_values
        valid_multidims[dim.getAttribute('name').lower()]['end'] = end_values
        valid_multidims[dim.getAttribute('name').lower()]['indexes'] = indexes

        if len(set([len(start_values), len(end_values), len(indexes)])) != 1:
            raise InternalError(
                'XML schema defines a multidimension keyword with start, end and indexes values of unequal length'
            )

    # Used to suggest valid keywords if an unknown one is found
    valid_invars_list = list(
        set([i.lower() for i in valid_dims] + [i.lower() for i in valid_multidims] + [i.lower() for i in valid_kws]) -
        set(blocked_kws)
    )

    # =================== Check for blocked keywords ===========================
    for keyword in input_params_internal:
        if keyword in blocked_kws:
            err_str = f"You should not provide explicitly keyword '{keyword}'."
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

    # from 5.0.2, this CANNOT be specified anymore!
    if Version(version) < Version('5.0.2'):
        # To be sure that things are read in angstrom - not possible in recent
        # versions
        input_params_internal['a'] = 1.

    # Get info on atomic species from the StructureData object
    atomic_species_list = [k.name for k in structure.kinds]

    try:
        calculation_type = input_params_internal['calculation']
    except KeyError as exception:
        raise QEInputValidationError(
            'Error, you need to specify at least the calculation type (among '
            f'{", ".join(list(valid_calculations_and_opt_namelists.keys()))})'
        ) from exception

    try:
        opt_namelists = valid_calculations_and_opt_namelists[calculation_type]
    except KeyError as exception:
        raise QEInputValidationError(
            f'Error, {calculation_type} is not a valid value for '
            f'the calculation type (valid values: {", ".join(list(valid_calculations_and_opt_namelists.keys()))}'
        ) from exception

    internal_dict = {i: {} for i in compulsory_namelists + opt_namelists}
    all_namelists = set(compulsory_namelists)
    for namelists in valid_calculations_and_opt_namelists.values():
        all_namelists.update(namelists)

    if not flat_mode:
        # Unexpected namelists specified by the user
        additional_namelists = sorted(all_input_namelists - set(all_namelists))
        if additional_namelists:
            formatted_namelists = ', '.join(additional_namelists)
            err_str = f'Error, the following namelists were specified but are not expected: {formatted_namelists}'
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

    # Empty list that contains the list of provided kws to check for
    # the compulsory ones at the end
    inserted_kws = []
    # I parse each element of the input dictionary
    for keyword, value in input_params_internal.items():
        keyword = keyword.lower()

        if keyword in valid_kws:
            # It is a variable
            found_var = valid_kws[keyword]
            namelist_name = found_var['namelist']
            if not flat_mode:
                input_namelist_name = input_original_namelists[keyword]
                if namelist_name != input_namelist_name:
                    err_str = (
                        f"Error, keyword '{keyword}' specified in namelist '{input_namelist_name}', but it should be "
                        f"instead in '{namelist_name}'"
                    )
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
            try:
                internal_dict[namelist_name][keyword] = _check_and_convert(keyword, value, found_var['expected_type'])
            except KeyError as exception:
                if namelist_name in all_namelists:
                    err_str = f'Error, namelist {namelist_name} not valid for calculation type {calculation_type}'
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str) from exception
                    else:
                        errors_list.append(err_str)
                else:
                    err_str = f'Error, unknown namelist {namelist_name}'
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str) from exception
                    else:
                        errors_list.append(err_str)
            except TypeError as exception:
                if stop_at_first_error:
                    raise
                else:
                    errors_list.append(str(exception))

        elif keyword in valid_dims:
            # It is an array
            found_var = valid_dims[keyword]
            namelist_name = found_var['namelist']
            if not flat_mode:
                input_namelist_name = input_original_namelists[keyword]
                if namelist_name != input_namelist_name:
                    err_str = (
                        f"Error, keyword '{keyword}' specified in namelist '{input_namelist_name}', but it should be "
                        f"instead in '{namelist_name}'"
                    )
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
            # I accept only ntyp or an integer as end_val
            if found_var['end_val'] == 'ntyp':
                if not isinstance(value, dict):
                    err_str = f"Error, expected dictionary to associate each specie to a value for keyword '{keyword}'."
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue

                outdict = {}
                for kindname, found_item in value.items():
                    if kindname not in atomic_species_list:
                        err_str = f"Error, '{kindname}' is not a valid kind name."
                        if stop_at_first_error:
                            raise QEInputValidationError(err_str)
                        else:
                            errors_list.append(err_str)
                            continue
                    try:
                        outdict[kindname] = _check_and_convert(keyword, found_item, found_var['expected_type'])
                    except TypeError as exception:
                        if stop_at_first_error:
                            raise
                        else:
                            errors_list.append(str(exception))

                try:
                    internal_dict[namelist_name][keyword] = outdict
                except KeyError as exception:
                    err_str = f'Error, unknown namelist {namelist_name}'
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str) from exception
                    else:
                        errors_list.append(err_str)
                        continue
            else:
                try:
                    end_value = int(found_var['end_val'])
                except ValueError as exception:
                    err_str = f"Error, invalid end value '{found_var['end_val']}' for keyword '{keyword}'."
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str) from exception
                    else:
                        errors_list.append(err_str)
                        continue
                if not isinstance(value, list) or len(value) != end_value:
                    err_str = f"Error, expecting a list of length {end_value} for keyword ' {keyword}'."
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue

                outlist = []
                for found_item in value:
                    if found_item is None:
                        # skip if the value is None (i.e., not provided)
                        outlist.append(None)
                    else:
                        try:
                            outlist.append(_check_and_convert(keyword, found_item, found_var['expected_type']))
                        except TypeError as exception:
                            if stop_at_first_error:
                                raise
                            else:
                                errors_list.append(str(exception))
                                outlist.append(None)

                try:
                    internal_dict[namelist_name][keyword] = outlist
                except KeyError as exception:
                    err_str = f'Error, unknown namelist {namelist_name}'
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str) from exception
                    else:
                        errors_list.append(err_str)
                        continue

        elif keyword in valid_multidims:
            # It is a multidimensional array

            variable = valid_multidims[keyword]
            indexes = variable['indexes']
            namelist_name = variable['namelist']

            # Create empty list for this keyword in the correct namelist
            try:
                internal_dict[namelist_name][keyword] = []
            except KeyError as exception:
                err_str = f'Error, unknown namelist {namelist_name}'
                if stop_at_first_error:
                    raise QEInputValidationError(err_str) from exception
                else:
                    errors_list.append(err_str)
                    continue

            for array in value:

                # Append empty list for this array of values
                internal_dict[namelist_name][keyword].append([])

                # Each array should contain N + 1 values, where N is the number of indexes for this multidimensional
                if len(array) != len(indexes) + 1:
                    err_str = f'Expected {len(indexes) + 1} values per array for kw {keyword}, got only {len(array)}.'
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue

                actual_value = array[-1]

                for i, index in enumerate(indexes):

                    index_value = array[i]

                    try:
                        int(variable['start'][i])
                    except ValueError as exception:
                        err_str = f"Error, invalid start value '{variable['start'][i]}' for keyword '{keyword}'."
                        if stop_at_first_error:
                            raise QEInputValidationError(err_str) from exception
                        else:
                            errors_list.append(err_str)
                            continue

                    end_value = variable['end'][i]

                    if end_value == 'ntyp':

                        kindname = index_value

                        if kindname not in atomic_species_list:
                            err_str = f"Error, '{kindname}' is not a valid kind name."
                            if stop_at_first_error:
                                raise QEInputValidationError(err_str)
                            else:
                                errors_list.append(err_str)
                                continue

                        internal_dict[namelist_name][keyword][-1].append(kindname)
                    else:
                        # Other types are assumed to be an integer
                        try:
                            index_value = int(index_value)
                        except ValueError as exception:
                            err_str = f'Error, only integer types are supported for index {index}, got {index_value}'
                            if stop_at_first_error:
                                raise QEInputValidationError(err_str) from exception
                            else:
                                errors_list.append(err_str)
                                continue

                        internal_dict[namelist_name][keyword][-1].append(index_value)

                # Validate the actual value, convert it and append to the current array
                converted = _check_and_convert(keyword, actual_value, variable['expected_type'])
                internal_dict[namelist_name][keyword][-1].append(converted)

        else:
            # Neither a variable nor an array
            err_str = f'Problem parsing keyword {keyword}. '
            similar_kws = difflib.get_close_matches(keyword, valid_invars_list)
            if len(similar_kws) == 1:
                err_str += f'Maybe you wanted to specify {similar_kws[0]}?'
            elif len(similar_kws) > 1:
                err_str += f"Maybe you wanted to specify one of these: {', '.join(similar_kws)}?"
            else:
                err_str += '(No similar keywords found...)'
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

        # Used to check if all compulsory variables are set
        inserted_kws += [keyword]

    # ============== I check here compulsory variables ===========
    missing_kws = compulsory_kws - set(inserted_kws)
    if missing_kws:
        err_str = f"Missing compulsory variables: {', '.join(missing_kws)}."
        if stop_at_first_error:
            raise QEInputValidationError(err_str)
        else:
            errors_list.append(err_str)

    if errors_list:
        raise QEInputValidationError(f'Errors! {len(errors_list)} issues found:\n* ' + '\n* '.join(errors_list))

    return internal_dict
