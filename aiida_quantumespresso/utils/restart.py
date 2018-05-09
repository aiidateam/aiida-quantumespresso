# -*- coding: utf-8 -*-
"""
Utility functions to return builders ready to be submitted 
for restarting a Quantum ESPRESSO calculation (or to apply small 
modifications).
"""
from aiida.common.datastructures import calc_states
import aiida_quantumespresso.calculations.ph
from aiida.common.exceptions import InputValidationError 
from aiida.common.links import LinkType

def create_restart_ph(parent_calc, force_restart=False,
                      parent_folder_symlink=aiida_quantumespresso.calculations.ph._default_symlink_usage):
    """
    This function creates a builder to restart a ph.x Quantum ESPRESSO 
    calculation that was not completed before (like max walltime reached...). This means that it might not be useful to restart a FAILED calculation.

    It returns a `builder` for a new calculation, with all links prepared 
    but not stored in DB.
    To submit it, simply call::

        from aiida.work.launch import submit
        submit(builder)

    :param bool force_restart: restart also if parent is not in FINISHED 
        state (e.g. FAILED, IMPORTED, etc.). Default=False.
    :param bool parent_folder_symlink: sets the value of the 
        `PARENT_FOLDER_SYMLINK` in the `settings` input.
    """
    from aiida.orm import DataFactory

    ParameterData = DataFactory('parameter')
    RemoteData = DataFactory('remote')

    if not isinstance(parent_calc, aiida_quantumespresso.calculations.ph.PhCalculation):
        raise TypeError(
            "This function can only deal with restarts of PhCalculations (QE ph.x codes), "
            "but I got {} instead".format(parent_calc.__class__.__name__))

    if parent_calc.get_state() != calc_states.FINISHED:
        if force_restart:
            pass
        else:
            raise InputValidationError(
                "Calculation to be restarted must be "
                "in the {} state. Otherwise, use the force_restart "
                "flag".format(calc_states.FINISHED) )
    
    inp = parent_calc.get_inputs_dict(link_type=LinkType.INPUT)
    code = inp['code']
    qpoints = inp['qpoints']
    
    old_inp_dict = inp['parameters'].get_dict()
    # add the restart flag
    old_inp_dict['INPUTPH']['recover'] = True
    inp_dict = ParameterData(dict=old_inp_dict) 

    remote_folders = parent_calc.get_outputs(type=RemoteData)
    if len(remote_folders)!=1:
        raise InputValidationError("More than one output RemoteData found "
                                    "in calculation {}".format(parent_calc.pk))
    remote_folder = remote_folders[0]
    
    builder = parent_calc.__class__.get_builder()
    
    if not 'Restart' in parent_calc.label:
        builder.label = (parent_calc.label + " Restart of {} {}.".format(
                                    parent_calc.__class__.__name__,parent_calc.pk)).lstrip()
    else:
        builder.label = ("Restart of {} {}.".format(parent_calc.__class__.__name__,parent_calc.pk)).lstrip()
    
    # set the parameters, and the (same) code and q-points
    builder.parameters = ParameterData(inp_dict)
    builder.code = code
    builder.qpoints = qpoints

    try:
        old_settings_dict = inp['settings'].get_dict()
    except KeyError:
        old_settings_dict = {}
    if parent_folder_symlink:
        old_settings_dict['PARENT_FOLDER_SYMLINK'] = True
        
    if old_settings_dict: # if not empty dictionary
        settings = ParameterData(dict=old_settings_dict)
        builder.settings = settings
    
    builder.parent_folder = remote_folder
    
    return builder