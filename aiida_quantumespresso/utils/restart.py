# -*- coding: utf-8 -*-
"""
Utility functions to return builders ready to be submitted 
for restarting a Quantum ESPRESSO calculation (or to apply small 
modifications).
"""
from __future__ import absolute_import
from aiida.common import InputValidationError 
from aiida.common.links import LinkType

def clone_calculation(calculation):
    """
    Create a clone of a Calculation node.
    This is a temporary workaround, it's better to create a new builder.

    :returns: an unstored clone of this node. Note: links are not copied
    """
    clone = self.__class__()
    clone.dbnode.dbcomputer = self._dbnode.dbcomputer
    clone.dbnode.type = self._dbnode.type
    clone.label = self.label
    clone.description = self.description

    for key, value in self.iterattrs():
        if key != Sealable.SEALED_KEY:
            clone.set_attribute(key, value)

    for path in self.get_folder_list():
        clone.add_path(self.get_abs_path(path), path)

    return clone

def _create_restart_pw_cp(parent_calc, force_restart, parent_folder_symlink,
                   use_output_structure,
                   restart_from_beginning):
    """
    Function to restart a pw.x or cp.x calculation that was not completed before
    (like max walltime reached...) i.e. not to restart a really FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared
    but not stored in DB.
    To submit it, simply call::

        from aiida.engine import submit
        submit(builder)

    :param parent_calc: the calculation you want to restart
    :param bool force_restart: restart also if parent is not in FINISHED
       state (e.g. FAILED, IMPORTED, etc.).
    :param bool parent_folder_symlink: if True, symlinks are used
       instead of hard copies of the files.
       Pass None for the default calculation behavior.
    :param bool use_output_structure: if True, the output structure
       of the restarted calculation is used if available, rather than its
       input structure. Useful for nscf or bands calculations, but it
       SHOULD NOT be used for the restart of a relaxation run.
    :param bool restart_from_beginning: If set to True, creates a copy of
       the parent calculation, without using the scratch for the restart.
       Useful to restart calculations that have crashed on the cluster for
       external reasons.
    """
    from aiida.plugins import DataFactory, CalculationFactory
    from aiida_quantumespresso.utils.pseudopotential import get_pseudos_of_calc

    Dict = DataFactory('dict')
    RemoteData = DataFactory('remote')

    if not isinstance(parent_calc, (CalculationFactory('quantumespresso.cp'), CalculationFactory('quantumespresso.pw'))):
        raise TypeError(
            "This function can only deal with restarts of PwCalculations or CpCalculations (QE pw.x or cp.x codes), "
            "but I got {} instead".format(parent_calc.__class__.__name__))

    # Check the calculation's state using ``from_attribute=True`` to
    # correctly handle IMPORTED calculations.
    if not parent_calc.is_finished_ok:
        if not force_restart:
            raise InputValidationError(
                "Calculation to be restarted must be finished ok. Otherwise, use the force_restart flag")

    inputs = parent_calc.get_incoming(link_type=LinkType.INPUT_CALC)

    if restart_from_beginning:
        inp_dict = inputs.get_node_by_label(parent_calc.get_linkname('parameters'))
    else:  # case of restart without using the parent scratch
        old_inp_dict = inputs.get_node_by_label(parent_calc.get_linkname('parameters')).get_dict()
        # add the restart flag
        old_inp_dict['CONTROL']['restart_mode'] = 'restart'
        inp_dict = Dict(dict=old_inp_dict)

    try:
        remote_folder = parent_calc.get_outcoming(node_class=RemoteData, link_type=LinkType.CREATE).one().node
    except ValueError:
        raise InputValidationError("No or more than one output RemoteData found in calculation {}".format(parent_calc.pk))

    builder = parent_calc.__class__.get_builder()

    # Set the same options
    # TODO: Implement this in a way, in AiiDA, to allow to get the values in a single call,
    #       returning either only those actually set, or all with their defaults,
    #       and here just call that method (to re-set those actually set).
    for option in builder.options._valid_fields:
        option_val = getattr(parent_calc, 'get_{}'.format(option))()
        if option_val:  # Skip None, empty strings, and empty dicts/lists
            setattr(builder.options, option, option_val)

    builder.label = parent_calc.label
    builder.description = "[Restart of {} {}]\n{}".format(
        parent_calc.__class__.__name__, parent_calc.uuid,
        parent_calc.description)

    # set the new links
    builder.parameters = inp_dict
    if use_output_structure:
        # use OUTPUT structure if available
        try:
            builder.structure = parent_calc.outputs.output_structure
        except AttributeError:
            builder.structure = inputs.get_node_by_label(parent_calc.get_linkname('structure'))
    else:
        builder.structure = inputs.get_node_by_label(parent_calc.get_linkname('structure'))

    # This distinguishes between pw.x and cp.x
    if parent_calc._use_kpoints:
        builder.kpoints = inputs.get_node_by_label(parent_calc.get_linkname('kpoints'))

    builder.code = inputs.get_node_by_label(parent_calc.get_linkname('code'))

    try:
        old_settings_dict = inputs.get_node_by_label(parent_calc.get_linkname('settings')).get_dict()
    except KeyError:
        old_settings_dict = {}
    if parent_folder_symlink is None:
        parent_folder_symlink = parent_calc._default_symlink_usage
    # Always set if it was already set. Otherwise, if it wasn't set, just set it if it's not the default
    if ('PARENT_FOLDER_SYMLINK' in old_settings_dict
            or parent_folder_symlink != parent_calc._default_symlink_usage):
        old_settings_dict['PARENT_FOLDER_SYMLINK'] = parent_folder_symlink

    if old_settings_dict:  # if not empty dictionary
        settings = Dict(dict=old_settings_dict)
        builder.settings = settings

    builder.parent_folder = remote_folder

    builder.pseudo = get_pseudos_of_calc(parent_calc)

    # Add also the vdw table, if the parent had one
    try:
        old_vdw_table = inputs.get_node_by_label(parent_calc.get_linkname('vdw_table'))
    except KeyError:
        # No VdW table
        pass
    else:
        builder.vdw_table = old_vdw_table

    return builder

def create_restart_pw(parent_calc, force_restart=False, parent_folder_symlink=None,
                   use_output_structure=False,
                   restart_from_beginning=False):
    """
    Function to restart a pw.x calculation that was not completed before
    (like max walltime reached...) i.e. not to restart a really FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared
    but not stored in DB.
    To submit it, simply call::

        from aiida.engine import submit
        submit(builder)

    :param parent_calc: the calculation you want to restart
    :param bool force_restart: restart also if parent is not in FINISHED
       state (e.g. FAILED, IMPORTED, etc.). Default=False.
    :param bool parent_folder_symlink: if True, symlinks are used
       instead of hard copies of the files. Default given by
       self._default_symlink_usage.
    :param bool use_output_structure: if True, the output structure
       of the restarted calculation is used if available, rather than its
       input structure. Useful for nscf or bands calculations, but it
       SHOULD NOT be used for the restart of a relaxation run.
       Default=False.
    :param bool restart_from_beginning: If set to True, creates a copy of
       the parent calculation, without using the scratch for the restart.
       Useful to restart calculations that have crashed on the cluster for
       external reasons. Default=False
    """
    return _create_restart_pw_cp(
        parent_calc=parent_calc, force_restart=force_restart,
        parent_folder_symlink=parent_folder_symlink,
        use_output_structure=use_output_structure,
        restart_from_beginning=restart_from_beginning)


def create_restart_cp(parent_calc, force_restart=False, parent_folder_symlink=None,
                   use_output_structure=False,
                   restart_from_beginning=False):
    """
    Function to restart a cp.x calculation that was not completed before
    (like max walltime reached...) i.e. not to restart a really FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared
    but not stored in DB.
    To submit it, simply call::

        from aiida.engine import submit
        submit(builder)

    :param parent_calc: the calculation you want to restart
    :param bool force_restart: restart also if parent is not in FINISHED
       state (e.g. FAILED, IMPORTED, etc.). Default=False.
    :param bool parent_folder_symlink: if True, symlinks are used
       instead of hard copies of the files. Default given by
       self._default_symlink_usage.
    :param bool use_output_structure: if True, the output structure
       of the restarted calculation is used if available, rather than its
       input structure. Useful for nscf or bands calculations, but it
       SHOULD NOT be used for the restart of a relaxation run.
       Default=False.
    :param bool restart_from_beginning: If set to True, creates a copy of
       the parent calculation, without using the scratch for the restart.
       Useful to restart calculations that have crashed on the cluster for
       external reasons. Default=False
    """
    return _create_restart_pw_cp(
        parent_calc=parent_calc, force_restart=force_restart,
        parent_folder_symlink=parent_folder_symlink,
        use_output_structure=use_output_structure,
        restart_from_beginning=restart_from_beginning)

def create_restart_neb(parent_calc, force_restart=False, parent_folder_symlink=None):
    """
    Function to restart a calculation that was not completed before
    (like max walltime reached...) i.e. not to restart a really FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared
    but not stored in DB.
    To submit it, simply call::

        from aiida.engine import submit
        submit(builder)

    :param parent_calc: the calculation you want to restart
    :param bool force_restart: restart also if parent is not in FINISHED
    state (e.g. FAILED, IMPORTED, etc.). Default=False.
    :param bool parent_folder_symlink: if True, symlinks are used
       instead of hard copies of the files.
       Pass None for the default calculation behavior.
    """
    from aiida.plugins import DataFactory, CalculationFactory
    from aiida_quantumespresso.utils.pseudopotential import get_pseudos_of_calc

    Dict = DataFactory('dict')
    RemoteData = DataFactory('remote')

    if not isinstance(parent_calc, CalculationFactory('quantumespresso.neb')):
        raise TypeError(
            "This function can only deal with restarts of NebCalculations (QE neb.x), "
            "but I got {} instead".format(parent_calc.__class__.__name__))

    # Check the calculation's state using ``from_attribute=True`` to
    # correctly handle IMPORTED calculations.
    if not force_restart and not parent_calc.is_finished_ok:
        raise exceptions.InputValidationError(
            "Calculation to be restarted must be finshed ok. Otherwise, use the force_restart flag")

    inputs = parent_calc.get_incoming(link_type=LinkType.INPUT_CALC)

    old_inp_dict = inputs.get_node_by_label(parent_calc.get_linkname('neb_parameters')).get_dict()
    # add the restart flag
    old_inp_dict['PATH']['restart_mode'] = 'restart'
    inp_dict = Dict(dict=old_inp_dict)

    try:
        remote_folder = parent_calc.get_outcoming(node_class=RemoteData, link_type=LinkType.CREATE).one().node
    except ValueError:
        raise InputValidationError("No or more than one output RemoteData found in calculation {}".format(parent_calc.pk))

    builder = parent_calc.__class__.get_builder()

    # Set the same options
    # TODO: Implement this in a way, in AiiDA, to allow to get the values in a single call,
    #       returning either only those actually set, or all with their defaults,
    #       and here just call that method (to re-set those actually set).
    for option in builder.options._valid_fields:
        option_val = getattr(parent_calc, 'get_{}'.format(option))()
        if option_val:  # Skip None, empty strings, and empty dicts/lists
            setattr(builder.options, option, option_val)

    builder.label = parent_calc.label
    builder.description = "[Restart of {} {}]\n{}".format(
        parent_calc.__class__.__name__, parent_calc.uuid,
        parent_calc.description)

    # set the new links
    builder.neb_parameters = inp_dict

    builder.pw_parameters = inputs.get_node_by_label(parent_calc.get_linkname('pw_parameters'))

    builder.first_structure = inputs.get_node_by_label(parent_calc.get_linkname('first_structure'))
    builder.last_structure = inputs.get_node_by_label(parent_calc.get_linkname('last_structure'))

    if parent_calc._use_kpoints:
        builder.kpoints = inputs.get_node_by_label(parent_calc.get_linkname('kpoints'))
    builder.code = inputs.get_node_by_label(parent_calc.get_linkname('code'))

    try:
        old_settings_dict = inputs.get_node_by_label(parent_calc.get_linkname('settings')).get_dict()
    except KeyError:
        old_settings_dict = {}
    if parent_folder_symlink is None:
        parent_folder_symlink = parent_calc._default_symlink_usage
    # Always set if it was already set. Otherwise, if it wasn't set, just set it if it's not the default
    if ('PARENT_FOLDER_SYMLINK' in old_settings_dict
            or parent_folder_symlink != parent_calc._default_symlink_usage):
        old_settings_dict['PARENT_FOLDER_SYMLINK'] = parent_folder_symlink

    if old_settings_dict:  # if not empty dictionary
        settings = Dict(dict=old_settings_dict)
        builder.settings = settings

    builder.parent_folder = remote_folder

    builder.pseudo = get_pseudos_of_calc(parent_calc)

    # Add also the vdw table, if the parent had one
    try:
        old_vdw_table = inputs.get_node_by_label(parent_calc.get_linkname('vdw_table'))
    except KeyError:
        # No VdW table
        pass
    else:
        builder.vdw_table = old_vdw_table

    return builder


def create_restart_ph(parent_calc, force_restart=False,
                      parent_folder_symlink=None):
    """
    This function creates a builder to restart a ph.x Quantum ESPRESSO 
    calculation that was not completed before (like max walltime reached...).
    This means that it might not be useful to restart a FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared 
    but not stored in DB.
    To submit it, simply call::

        from aiida.engine import submit
        submit(builder)

    :param parent_calc: the calculation you want to restart
    :param bool force_restart: restart also if parent is not in FINISHED 
        state (e.g. FAILED, IMPORTED, etc.). Default=False.
    :param bool parent_folder_symlink: sets the value of the 
        `PARENT_FOLDER_SYMLINK` in the `settings` input. Default: the value specified
        in aiida_quantumespresso.calculations.ph._default_symlink_usage.
    """
    from aiida.plugins import DataFactory, CalculationFactory

    Dict = DataFactory('dict')
    RemoteData = DataFactory('remote')

    if not isinstance(parent_calc, CalculationFactory('quantumespresso.ph')):
        raise TypeError(
            "This function can only deal with restarts of PhCalculations (QE ph.x codes), "
            "but I got {} instead".format(parent_calc.__class__.__name__))

    if not force_restart and not parent_calc.is_finished_ok:
        raise exceptions.InputValidationError(
            "Calculation to be restarted must be finshed ok. Otherwise, use the force_restart flag")
    
    inputs = parent_calc.get_incoming(link_type=LinkType.INPUT_CALC)
    code = inputs.get_node_by_label('code')
    qpoints = inputs.get_node_by_label('qpoints')
    
    old_inp_dict = inputs.get_node_by_label('parameters').get_dict()
    # add the restart flag
    old_inp_dict['INPUTPH']['recover'] = True
    inp_dict = Dict(dict=old_inp_dict) 

    try:
        remote_folder = parent_calc.get_outcoming(node_class=RemoteData, link_type=LinkType.CREATE).one().node
    except ValueError:
        raise InputValidationError("No or more than one output RemoteData found in calculation {}".format(parent_calc.pk))
    
    builder = parent_calc.__class__.get_builder()

    # Set the same options
    # TODO: Implement this in a way, in AiiDA, to allow to get the values in a single call,
    #       returning either only those actually set, or all with their defaults,
    #       and here just call that method (to re-set those actually set).
    for option in builder.options._valid_fields:
        option_val = getattr(parent_calc, 'get_{}'.format(option))()
        if option_val:  # Skip None, empty strings, and empty dicts/lists
            setattr(builder.options, option, option_val)

    builder.label = parent_calc.label
    builder.description = "[Restart of {} {}]\n{}".format(
        parent_calc.__class__.__name__, parent_calc.uuid,
        parent_calc.description)
    
    # set the parameters, and the (same) code and q-points
    builder.parameters = Dict(inp_dict)
    builder.code = code
    builder.qpoints = qpoints

    try:
        old_settings_dict = inputs.get_node_by_label('settings').get_dict()
    except KeyError:
        old_settings_dict = {}
    if parent_folder_symlink is None:
        parent_folder_symlink = parent_calc._default_symlink_usage
    # Always set if it was already set. Otherwise, if it wasn't set, just set it if it's not the default
    if ('PARENT_FOLDER_SYMLINK' in old_settings_dict
            or parent_folder_symlink != parent_calc._default_symlink_usage):
        old_settings_dict['PARENT_FOLDER_SYMLINK'] = parent_folder_symlink
        
    if old_settings_dict: # if not empty dictionary
        settings = Dict(dict=old_settings_dict)
        builder.settings = settings
    
    builder.parent_folder = remote_folder
    
    return builder
