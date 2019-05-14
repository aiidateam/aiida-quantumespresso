# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBandsWorkChain` for testing and demonstration purposes."""
from __future__ import absolute_import
import click

from aiida.cmdline.params import options as options_cli
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from ...cli import workflow_launch
from ...utils import launch
from ...utils import options as options_qe
from ...utils import validate


@workflow_launch.command('pw-bands')
@options_cli.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(required=True)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_DISTANCE()
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.HUBBARD_U()
@options_qe.HUBBARD_V()
@options_qe.HUBBARD_FILE()
@options_qe.STARTING_MAGNETIZATION()
@options_qe.SMEARING()
@options_qe.AUTOMATIC_PARALLELIZATION()
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(code, structure, pseudo_family, kpoints_distance, ecutwfc, ecutrho, hubbard_u, hubbard_v,
                    hubbard_file_pk, starting_magnetization, smearing, automatic_parallelization, clean_workdir,
                    max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a `PwBandsWorkChain`."""
    from aiida.orm import Bool, Float, Str, Dict
    from aiida.plugins import WorkflowFactory
    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    parameters = {
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
        },
    }

    try:
        hubbard_file = validate.validate_hubbard_parameters(structure, parameters, hubbard_u, hubbard_v,
                                                            hubbard_file_pk)
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    try:
        validate.validate_starting_magnetization(structure, parameters, starting_magnetization)
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    try:
        validate.validate_smearing(parameters, smearing)
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    pseudo_family = Str(pseudo_family)
    parameters = Dict(dict=parameters)

    inputs = {
        'structure': structure,
        'relax': {
            'base': {
                'code': code,
                'pseudo_family': pseudo_family,
                'kpoints_distance': Float(kpoints_distance),
                'parameters': parameters,
            },
            'meta_convergence': Bool(False),
        },
        'scf': {
            'code': code,
            'pseudo_family': pseudo_family,
            'kpoints_distance': Float(kpoints_distance),
            'parameters': parameters,
        },
        'bands': {
            'code': code,
            'pseudo_family': pseudo_family,
            'kpoints_distance': Float(kpoints_distance),
            'parameters': parameters,
        }
    }

    if hubbard_file:
        inputs['relax']['base']['hubbard_file'] = hubbard_file
        inputs['scf']['base']['hubbard_file'] = hubbard_file
        inputs['bands']['base']['hubbard_file'] = hubbard_file

    if automatic_parallelization:
        auto_para = Dict(dict=get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds))
        inputs['relax']['base']['automatic_parallelization'] = auto_para
        inputs['scf']['automatic_parallelization'] = auto_para
        inputs['bands']['automatic_parallelization'] = auto_para
    else:
        options = Dict(dict=get_default_options(max_num_machines, max_wallclock_seconds, with_mpi))
        inputs['relax']['base']['options'] = options
        inputs['scf']['options'] = options
        inputs['bands']['options'] = options

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    launch.launch_process(WorkflowFactory('quantumespresso.pw.bands'), daemon, **inputs)
