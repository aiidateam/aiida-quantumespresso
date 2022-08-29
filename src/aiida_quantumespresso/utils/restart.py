# -*- coding: utf-8 -*-
"""Utility functions to return process builders ready to be submitted for restarting a Quantum ESPRESSO calculation."""
from aiida_quantumespresso.calculations.cp import CpCalculation
from aiida_quantumespresso.calculations.neb import NebCalculation
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation


def get_builder_restart(node, from_scratch=False, use_symlink=False):
    """Create a restart `ProcessBuilder` from a completed `CalcJobNode`.

    The restart builder will be ready to be launched but of course one can further update the inputs before doing so.
    To launch the restart calculation, simply run or submit it like you would normally::

        from aiida.engine import submit
        builder = get_builder_restart(node)
        submit(builder)

    :param node: the `CalcJobNode` instance from which to create a restart builder
    :param from_scratch: boolean, if True, will restart from scratch
    :param use_symlink: boolean, if True, will symlink the parent folder instead of copying the contents
    :return: a `ProcessBuilder` instance configured for a restart
    """
    from aiida.orm import Dict

    supported = (CpCalculation, NebCalculation, PhCalculation, PwCalculation)

    if node.process_class not in supported:
        raise TypeError(f'calculation class `{node.process_class}` of {node} is not one of {supported}')

    builder = node.get_builder_restart()

    if not from_scratch:
        builder.parent_folder = node.base.links.get_outgoing(link_label_filter='remote_folder').one().node
    else:
        builder.pop('parent_folder', None)

    # Update the parameters to set the correct restart flag
    parameters = builder.parameters.get_dict()

    if node.process_class in (CpCalculation, PwCalculation):
        parameters.setdefault('CONTROL', {})['restart_mode'] = 'from_scratch' if from_scratch else 'restart'
    elif node.process_class is NebCalculation:
        parameters.setdefault('PATH', {})['restart_mode'] = 'from_scratch' if from_scratch else 'restart'
    elif node.process_class is PhCalculation:
        parameters.setdefault('INPUTPH', {})['recover'] = True

    # If it was not already set, use the value passed as an argument or fallback to the class default
    settings = builder.settings.get_dict()
    settings.setdefault('PARENT_FOLDER_SYMLINK', use_symlink or node.process_class._default_symlink_usage)  # pylint: disable=protected-access

    builder.parameters = Dict(parameters)
    builder.settings = Dict(settings)

    return builder
