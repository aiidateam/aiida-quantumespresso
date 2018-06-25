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
    help='the plane wave cutoff energy in Ry'
)

ecutrho = overridable_option(
    '-R', '--ecutrho', type=click.FLOAT, default=240., show_default=True,
    help='the charge density cutoff energy in Ry'
)

hubbard_u = overridable_option(
    '-U', '--hubbard-u', nargs=2, multiple=True, type=click.Tuple([unicode, float]),
    help='add a Hubbard U term to a specific kind', metavar='<KIND MAGNITUDE>...'
)

hubbard_v = overridable_option(
    '-V', '--hubbard-v', nargs=4, multiple=True, type=click.Tuple([int, int, int, float]),
    help='add a Hubbard V interaction between two sites', metavar='<SITE SITE TYPE MAGNITUDE>...'
)

hubbard_file = overridable_option(
    '-H', '--hubbard-file', 'hubbard_file_pk', type=click.INT,
    help='the pk of a SinglefileData containing Hubbard parameters from a HpCalculation to use as input for Hubbard V'
)

starting_magnetization = overridable_option(
    '-M', '--starting-magnetization', nargs=2, multiple=True, type=click.Tuple([unicode, float]),
    help='add a starting magnetization to a specific kind', metavar='<KIND MAGNITUDE>...'
)

smearing = overridable_option(
    '-S', '--smearing', nargs=2, default=(None, None), type=click.Tuple([unicode, float]),
    help='add smeared occupations by specifying the type and amount of smearing', metavar='<TYPE DEGAUSS>'
)
