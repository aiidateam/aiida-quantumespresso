# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import xml.dom.minidom
import os
import difflib
import copy
from aiida.common import InputValidationError, InternalError
# Can also try to use LooseVersion instead, if more complicated things are
# required, e.g. with strings. But be careful, check if the behavior in
# this case is the intended one.
from distutils.version import StrictVersion
import six

class QEInputValidationError(InputValidationError):
    """
    This class is the exception that is generated by the parser when it
    encounters an error while creating the input file of Quantum ESPRESSO.
    """
    pass


def _check_and_convert(kw, val, expected_type):
    """
    val: the value to be read and converted to a Fortran-friendly string.
    expected_type: a string with the expected type. Can be:
      INTEGER
      REAL
      CHARACTER
      LOGICAL
    """

    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    if expected_type.upper() == 'LOGICAL':
        if isinstance(val, bool):
            outval = val
        else:
            raise TypeError(
                'Expected a boolean for keyword {}, found {} instead'.format(
                kw, type(val)))
    elif expected_type.upper() == 'REAL':
        if isinstance(val, six.integer_types):
            outval = float(val)
        elif isinstance(val, float):
            outval = val
        else:
            raise TypeError(
                'Expected a float for keyword {}, found {} instead'.format(
                kw, type(val)))
    elif expected_type.upper() == 'INTEGER':
        if isinstance(val, six.integer_types):
            outval = val
        else:
            raise TypeError(
                'Expected an integer for keyword {}, found {} instead'.format(
                kw, type(val)))
    elif expected_type.upper() == 'CHARACTER':
        if isinstance(val, six.string_types):
            outval = val
        else:
            raise TypeError(
                'Expected a string for keyword {}, found {} instead'.format(
                kw, type(val)))
    else:
        raise InternalError('Unexpected type check for keyword {}: {})'.format(
            kw, expected_type.upper()))

    return outval

def pw_input_helper(input_params, structure,
    stop_at_first_error=False, flat_mode=False, version='6.2'):
    """
    Validate if the input dictionary for Quantum ESPRESSO is valid.
    Return the dictionary (possibly with small variations: e.g. convert
    integer to float where necessary, recreate the proper structure
    if flat_mode is True, ...) to use as input parameters (use_parameters)
    for the pw.x calculation.

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
    errors_list = []

    # =========== LIST OF KNOWN NAMELISTS, CARDS, VARIABLES, ... ===============
    compulsory_namelists = ['CONTROL', 'SYSTEM', 'ELECTRONS']

    valid_calculations_and_opt_namelists = {
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
        for nl, content in six.iteritems(input_params):
            if not isinstance(content, dict):
                raise QEInputValidationError(
                    "The content associated to the namelist '{}' must be a "
                    'dictionary'.format(nl))
            all_input_namelists.add(nl)
            for k, v in six.iteritems(content):
                input_params_internal[k] = copy.deepcopy(v)
                if k in input_original_namelists:
                    err_str = "The keyword '{}' was specified both in the "
                    'namelist {} and {}.'.format(k,
                        input_original_namelists[k], nl)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                input_original_namelists[k] = nl

    # List of the keywords that must not appear in the input
    # (e.g. because they are automatically filled in by the plugin)
    blocked_kws = [i.lower() for i in
                   ['pseudo_dir',
                    'outdir',
                    'ibrav',
                    'celldm',
                    'nat',
                    'ntyp',
                    'prefix',
                    'a', 'b', 'c', 'cosab', 'cosac', 'cosbc',
                    ]
                   ]
    # TODO: possibly add here above also restart_mode?

    # List of the keywords that must ALWAYS appear in the input
    compulsory_kws = {i.lower() for i in
                   ['calculation',
                    'ecutwfc',
                     ]
                   }

    # ===================== PARSING OF THE XML DEFINITION FILE ===============
    module_dir = os.path.dirname(__file__)
    if module_dir == '':
        module_dir = os.curdir
    xml_path = os.path.join(module_dir, 'INPUT_PW-{}.xml'.format(version))
    try:
        with open(xml_path, 'r') as f:
            dom = xml.dom.minidom.parse(f)
    except IOError:
        prefix = 'INPUT_PW-'
        suffix = '.xml'
        versions = [fname[len(prefix):-len(suffix)] for fname
                    in os.listdir(module_dir) if fname.startswith(prefix)
                    and fname.endswith(suffix)]
        versions = sorted(versions, key=lambda x: StrictVersion(x))
        strictversions = versions + [version]
        strictversions = sorted(strictversions, key=lambda x: StrictVersion(x))
        pos = strictversions.index(version)
        if pos == 0:
            add_str = ' (the version you specified is too old)'
        else:
            add_str = ' (the older, closest version you can use is {})'.format(
                strictversions[pos-1])
        raise QEInputValidationError(
            'Unknown Quantum Espresso version: {}. '
            'Available versions: {};{}'.format(version, ', '.join(versions),
            add_str))


    # ========== List of known PW variables (from XML file) ===============
    known_kws = dom.getElementsByTagName('var')
    valid_kws = {}
    for kw in known_kws:
        if kw in valid_kws:
            raise InternalError('Something strange, I found more than one '
                                "keyword '{}' in the XML description...".format(
                                kw))

        valid_kws[kw.getAttribute('name').lower()] = {}
        parent = kw
        try:
            while True:
                parent = parent.parentNode
                if parent.tagName == 'namelist':
                    valid_kws[kw.getAttribute('name').lower()]['namelist'] = \
                        parent.getAttribute('name').upper()
                    break
        except AttributeError:
            # There are also variables in cards instead of namelists:
            # I ignore them
            pass
                # raise QEInputValidationError("Unable to find namelist for "
                #     "keyword %s." % kw.getAttribute('name'))
        expected_type = kw.getAttribute('type')
        # Fix for groups of variables
        if expected_type == '':
            if kw.parentNode.tagName == 'vargroup':
                expected_type = kw.parentNode.getAttribute('type')
        valid_kws[kw.getAttribute('name').lower()]['expected_type'] = \
            expected_type.upper()


    # ====== List of known PW 'dimensions' (arrays) (from XML file) ===========
    known_dims = dom.getElementsByTagName('dimension')
    valid_dims = {}
    for dim in known_dims:
        if dim in valid_dims:
            raise InternalError('Something strange, I found more than one '
                "keyword '{}' in the XML description...".format(dim))

        valid_dims[dim.getAttribute('name').lower()] = {}
        parent = dim
        try:
            while True:
                parent = parent.parentNode
                if parent.tagName == 'namelist':
                    valid_dims[dim.getAttribute('name').lower()]['namelist'] = \
                        parent.getAttribute('name').upper()
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
        valid_dims[dim.getAttribute('name').lower()]['expected_type'] = \
            expected_type.upper()
        # I assume start_val is always 1
        start_val = dim.getAttribute('start')
        if start_val != '1':
            raise InternalError(
                "Wrong start value '{}' in input array (dimension) {}".format(
                    (start_val, dim.getAttribute('name'))))
        # I save the string as it is; somewhere else I will check for its value
        valid_dims[dim.getAttribute('name').lower()]['end_val'] = \
            dim.getAttribute('end')

    # Used to suggest valid keywords if an unknown one is found
    valid_invars_list = list(
        set([i.lower() for i in valid_dims.keys()] +
            [i.lower() for i in valid_kws.keys()]) - set(blocked_kws))

    # =================== Check for blocked keywords ===========================
    for kw in input_params_internal.keys():
        if kw in blocked_kws:
            err_str = "You should not provide explicitly keyword '{}'.".format(
                kw)
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

    # from 5.0.2, this CANNOT be specified anymore!
    if StrictVersion(version) < StrictVersion('5.0.2'):
        # To be sure that things are read in angstrom - not possible in recent
        # versions
        input_params_internal['a'] = 1.

    # Get info on atomic species from the StructureData object
    atomic_species_list = [k.name for k in structure.kinds]

    try:
        calculation_type = input_params_internal['calculation']
    except KeyError:
        raise QEInputValidationError('Error, you need to specify at least the '
            'calculation type (among {})'.format(
            ', '.join(list(valid_calculations_and_opt_namelists.keys()))))

    try:
        opt_namelists = valid_calculations_and_opt_namelists[calculation_type]
    except KeyError:
        raise QEInputValidationError('Error, {} is not a valid value for '
            'the calculation type (valid values: {})'.format(calculation_type,
            ', '.join(list(valid_calculations_and_opt_namelists.keys()))))

    internal_dict = {i: {} for i in compulsory_namelists + opt_namelists}
    all_namelists = set(compulsory_namelists)
    for namelists in valid_calculations_and_opt_namelists.values():
        all_namelists.update(namelists)

    if not flat_mode:
        # Unexpected namelists specified by the user
        additional_namelists = sorted(all_input_namelists - set(all_namelists))
        if additional_namelists:
            err_str = 'Error, the following namelists were specified but are ' \
                'not expected: {}'.format(
                    ', '.join(additional_namelists))
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

    # Empty list that contains the list of provided kws to check for
    # the compulsory ones at the end
    inserted_kws = []
    # I parse each element of the input dictionary
    for kw, value in six.iteritems(input_params_internal):
        #print kw, valid_kws[kw.lower()]
        kw = kw.lower()

        if kw in valid_kws:
            # It is a variable
            found_var = valid_kws[kw]
            namelist_name = found_var['namelist']
            if not flat_mode:
                input_namelist_name = input_original_namelists[kw]
                if namelist_name != input_namelist_name:
                    err_str = \
                        "Error, keyword '{}' specified in namelist '{}', " \
                        "but it should be instead in '{}'".format(
                            kw, input_namelist_name, namelist_name)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
            try:
                internal_dict[namelist_name][kw] = _check_and_convert(
                    kw, value, found_var['expected_type'])
            except KeyError:
                if namelist_name in all_namelists:
                    err_str = \
                        'Error, namelist {} not valid for calculation type ' \
                        '{}'.format(namelist_name, calculation_type)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                else:
                    err_str = 'Error, unknown namelist ' \
                        '{}'.format(namelist_name)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
            except TypeError as e:
                    if stop_at_first_error:
                        raise
                    else:
                        errors_list.append(e.message)

        elif kw in valid_dims:
            # It is an array
            found_var = valid_dims[kw]
            namelist_name = found_var['namelist']
            if not flat_mode:
                input_namelist_name = input_original_namelists[kw]
                if namelist_name != input_namelist_name:
                    err_str = \
                        "Error, keyword '{}' specified in namelist '{}', " \
                        "but it should be instead in '{}'".format(
                            kw, input_namelist_name, namelist_name)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
            ## I accept only ntyp or an integer as end_val
            if found_var['end_val'] == 'ntyp':
                if not isinstance(value, dict):
                    err_str = \
                        'Error, expecting a dictionary to associate each ' \
                        "specie to a value for keyword '{}'.".format(kw)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue

                outdict = {}
                for kindname, found_item in six.iteritems(value):
                    if kindname not in atomic_species_list:
                        err_str = \
                        "Error, '{}' is not a valid kind name.".format(kindname)
                        if stop_at_first_error:
                            raise QEInputValidationError(err_str)
                        else:
                            errors_list.append(err_str)
                            continue
                    try:
                        outdict[kindname] = _check_and_convert(kw, found_item,
                            found_var['expected_type'])
                    except TypeError:
                        if stop_at_first_error:
                            raise
                        else:
                            errors_list.append(e.message)

                try:
                    internal_dict[namelist_name][kw] = outdict
                except KeyError:
                    err_str = \
                        'Error, unknown namelist {}'.format(namelist_name)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue
            else:
                try:
                    end_value = int(found_var['end_val'])
                except ValueError:
                    err_str = \
                        "Error, invalid end value '{}' for keyword '{}'.".format(
                        (found_var['end_val'], kw))
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue
                if not isinstance(value, list) or len(value) != end_value:
                    err_str = \
                        'Error, expecting a list of length {} for keyword ' \
                        "'{}'.".format(end_value, kw)
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
                            outlist.append(_check_and_convert(kw, found_item,
                                found_var['expected_type']))
                        except TypeError as e:
                            if stop_at_first_error:
                                raise
                            else:
                                errors_list.append(e.message)
                                outlist.append(None)

                try:
                    internal_dict[namelist_name][kw] = outlist
                except KeyError:
                    err_str = \
                        'Error, unknown namelist {}'.format(namelist_name)
                    if stop_at_first_error:
                        raise QEInputValidationError(err_str)
                    else:
                        errors_list.append(err_str)
                        continue
        else:
            # Neither a variable nor an array
            err_str = 'Problem parsing keyword {}. '.format(kw)
            similar_kws = difflib.get_close_matches(kw, valid_invars_list)
            if len(similar_kws)==1:
                err_str += 'Maybe you wanted to specify {}?'.format(
                    similar_kws[0])
            elif len(similar_kws) > 1:
                err_str += 'Maybe you wanted to specify one of these: ' \
                    '{}?'.format(', '.join(similar_kws))
            else:
                err_str += '(No similar keywords found...)'
            if stop_at_first_error:
                raise QEInputValidationError(err_str)
            else:
                errors_list.append(err_str)

        # Used to check if all compulsory variables are set
        inserted_kws += [kw]

    # ============== I check here compulsory variables ===========
    missing_kws = compulsory_kws - set(inserted_kws)
    if len(missing_kws) != 0:
        err_str = 'Missing compulsory variables: {}.'.format(
            ', '.join(missing_kws))
        if stop_at_first_error:
            raise QEInputValidationError(err_str)
        else:
            errors_list.append(err_str)

    if errors_list:
        raise QEInputValidationError(
            'Errors! {} issues found:\n* '.format(len(errors_list)) +
            '\n* '.join(errors_list))

    return internal_dict


if __name__ == '__main__':
     # An example of usage
     from aiida.orm import load_node
     structure = DataFactory('structure')(cell=[[1,0,0],[0,1,0],[0,0,1]])
     structure.append_atom(symbols='Si', position=[0,0,0])
     structure.append_atom(symbols='O', position=[0.5,0.5,0.5])

     try:
         print(validate_pw_input({
             'calculation': 'vc-relax',
             'ecutwfc': 30.,
             'lda_plus_u': True,
             'lda_plus_u_kind': 2,
             'ion_temperature': 'a',
 #            'hubbard_u': [1, None],
             'hubbard_u': {'O': 1},
             },
             structure, flat_mode = True,
             version = '5.1'))
     except QEInputValidationError as e:
         print('*'*72)
         print('* ERROR !')
         print('*'*72)
         print(e.message)

     try:
         print(validate_pw_input(
             {
                 'CONTROL': {
                     'calculation': 'vc-relax'
                     },
                 'IONS': {
                      'ion_temperature': 'a'
                     },
                 'CELL': {
                     },
                 'ELECTRONS': {
                     },
                 'SYSTEM': {
                     'lda_plus_u_kind': 2,
                     'ecutwfc': 30.0,
                     'hubbard_u': {'O': 1.0},
                     'lda_plus_u': True}
             },
             structure, flat_mode = False))
     except QEInputValidationError as e:
         print('*'*72)
         print('* ERROR !')
         print('*'*72)
         print(e.message)


