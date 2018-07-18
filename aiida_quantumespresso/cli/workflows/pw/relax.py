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
@click.option('-k', '--kpoints-distance', type=click.FLOAT, required=True, default=0.15)
@click.option(
    '-f', '--final-scf', is_flag=True, default=False, show_default=True,
    help='run a final scf calculation for the final relaxed structure'
)
@click.option(
    '-g', '--group', type=click.STRING, required=False,
    help='the label of a Group to add the final PwCalculation to in case of success'
)
def launch(
    code, structure, pseudo_family, kpoints, max_num_machines, max_wallclock_seconds, daemon, ecutwfc, ecutrho,
    hubbard_u, hubbard_v, hubbard_file_pk, starting_magnetization, smearing, automatic_parallelization, clean_workdir,
    final_scf, group, kpoints_distance):
    """
    Run the PwRelaxWorkChain for a given input structure
    """
    from aiida.orm.data.base import Bool, Float, Str
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import WorkflowFactory
    from aiida.work.launch import run, submit
    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')

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

    inputs = {
        'structure': structure,
        'base': {
            'code': code,
            'pseudo_family': Str(pseudo_family),
            'kpoints_distance': Float(kpoints_distance),
            'parameters': ParameterData(dict=parameters),
        }
    }

    if automatic_parallelization:
        automatic_parallelization = get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds)
        inputs['base']['automatic_parallelization'] = ParameterData(dict=automatic_parallelization)
    else:
        options = get_default_options(max_num_machines, max_wallclock_seconds)
        inputs['base']['options'] = ParameterData(dict=options)

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if final_scf:
        inputs['final_scf'] = Bool(True)

    if group:
        inputs['group'] = Str(group)

    if daemon:
        workchain = submit(PwRelaxWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwRelaxWorkChain.__name__, workchain.pk))
    else:
        run(PwRelaxWorkChain, **inputs)
