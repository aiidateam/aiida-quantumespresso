# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options
from aiida_quantumespresso.utils.cli import options as options_qe
from aiida_quantumespresso.utils.cli import validate


@command()
@options.code(callback_kwargs={'entry_point': 'quantumespresso.pw'})
@options.structure()
@options.pseudo_family()
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
@options.daemon()
@options_qe.ecutwfc()
@options_qe.ecutrho()
@options_qe.hubbard_u()
@options_qe.hubbard_v()
@options_qe.hubbard_file()
@options_qe.starting_magnetization()
@options_qe.smearing()
@options_qe.automatic_parallelization()
@options_qe.clean_workdir()
def launch(
    code, structure, pseudo_family, kpoints, max_num_machines, max_wallclock_seconds, daemon, ecutwfc, ecutrho,
    hubbard_u, hubbard_v, hubbard_file_pk, starting_magnetization, smearing, automatic_parallelization, clean_workdir):
    """
    Run the PwBandsWorkChain for a given input structure
    """
    from aiida.orm.data.base import Bool, Float, Str
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import WorkflowFactory
    from aiida.work.launch import run, submit
    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    PwBandsWorkChain = WorkflowFactory('quantumespresso.pw.bands')

    parameters = {
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
        },
    }

    try:
        hubbard_file = validate.validate_hubbard_parameters(structure, parameters, hubbard_u, hubbard_v, hubbard_file_pk)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    try:
        validate.validate_starting_magnetization(structure, parameters, starting_magnetization)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    try:
        validate.validate_smearing(parameters, smearing)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    pseudo_family = Str(pseudo_family)
    parameters = ParameterData(dict=parameters)

    inputs = {
        'structure': structure,
        'relax': {
            'base': {
                'code': code,
                'pseudo_family': pseudo_family,
                'kpoints_distance': Float(1.0),
                'parameters': parameters,
                'meta_convergence': Bool(False),
            }
        },
        'scf': {
            'code': code,
            'pseudo_family': pseudo_family,
            'kpoints_distance': Float(1.0),
            'parameters': parameters,
        },
        'bands': {
            'code': code,
            'pseudo_family': pseudo_family,
            'kpoints_distance': Float(1.0),
            'parameters': parameters,
        }
    }

    if automatic_parallelization:
        auto_para = ParameterData(dict=get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds))
        inputs['relax']['base']['automatic_parallelization'] = auto_para
        inputs['scf']['automatic_parallelization'] = auto_para
        inputs['bands']['automatic_parallelization'] = auto_para
    else:
        options = ParameterData(dict=get_default_options(max_num_machines, max_wallclock_seconds))
        inputs['relax']['base']['options'] = options
        inputs['scf']['options'] = options
        inputs['bands']['options'] = options

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = submit(PwBandsWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandsWorkChain.__name__, workchain.pk))
    else:
        run(PwBandsWorkChain, **inputs)
