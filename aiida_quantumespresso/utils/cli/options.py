# -*- coding: utf-8 -*-
import click
from aiida.utils.cli.options import overridable_option


automatic_parallelization = overridable_option(
    '-a', '--automatic-parallelization', is_flag=True, default=False, show_default=True,
    help='enable the automatic parallelization option of the workchain'
)

clean_workdir = overridable_option(
    '-x', '--clean-workdir', is_flag=True, default=False, show_default=True,
    help='clean the remote folder of all the launched calculations after completion of the workchain'
)

ecutwfc = overridable_option(
    '-W', '--ecutwfc', type=click.FLOAT, default=30., show_default=True,
    help='The plane wave cutoff energy in Ry'
)

ecutrho = overridable_option(
    '-R', '--ecutrho', type=click.FLOAT, default=240., show_default=True,
    help='The charge density cutoff energy in Ry'
)

hubbard_u = overridable_option(
    '-U', '--hubbard-u', nargs=2, multiple=True, type=click.Tuple([unicode, float]),
    help='Add a Hubbard U term to a specific kind'
)