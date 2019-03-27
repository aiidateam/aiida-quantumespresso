# -*- coding: utf-8 -*-
"""Pre-defined overridable options for commonly used command line interface parameters."""
from __future__ import absolute_import

import six
import click

from aiida.cmdline.params import types
from aiida.cmdline.params.options import OverridableOption

from . import validate

STRUCTURE = OverridableOption(
    '-s', '--structure', type=types.DataParamType(sub_classes=('aiida.data:structure',)), help=u'StructureData node.')

PSEUDO_FAMILY = OverridableOption(
    '-p', '--pseudo-family', 'pseudo_family', type=click.STRING, help=u'Pseudo potential family name.')

KPOINTS_DISTANCE = OverridableOption(
    '-K',
    '--kpoints-distance',
    type=click.FLOAT,
    default=0.5,
    show_default=True,
    help=u'The minimal distance between k-points in reciprocal space in inverse Ångström.')

KPOINTS_MESH = OverridableOption(
    '-k',
    '--kpoints-mesh',
    'kpoints_mesh',
    nargs=3,
    type=click.INT,
    show_default=True,
    callback=validate.validate_kpoints_mesh,
    help=u'The number of points in the kpoint mesh along each basis vector.')

MAX_NUM_MACHINES = OverridableOption(
    '-m',
    '--max-num-machines',
    type=click.INT,
    default=1,
    show_default=True,
    help=u'The maximum number of machines (nodes) to use for the calculations.')

MAX_WALLCLOCK_SECONDS = OverridableOption(
    '-w',
    '--max-wallclock-seconds',
    type=click.INT,
    default=1800,
    show_default=True,
    help=u'the maximum wallclock time in seconds to set for the calculations.')

WITH_MPI = OverridableOption(
    '-i', '--with-mpi', is_flag=True, default=False, show_default=True, help=u'Run the calculations with MPI enabled.')

DAEMON = OverridableOption(
    '-d',
    '--daemon',
    is_flag=True,
    default=False,
    show_default=True,
    help=u'Submit the process to the daemon instead of running it locally.')

AUTOMATIC_PARALLELIZATION = OverridableOption(
    '-a',
    '--automatic-parallelization',
    is_flag=True,
    default=False,
    show_default=True,
    help=u'Enable the automatic parallelization option of the workchain.')

CLEAN_WORKDIR = OverridableOption(
    '-x',
    '--clean-workdir',
    is_flag=True,
    default=False,
    show_default=True,
    help=u'Clean the remote folder of all the launched calculations after completion of the workchain.')

ECUTWFC = OverridableOption(
    '-W', '--ecutwfc', type=click.FLOAT, default=30., show_default=True, help=u'The plane wave cutoff energy in Ry.')

ECUTRHO = OverridableOption(
    '-R',
    '--ecutrho',
    type=click.FLOAT,
    default=240.,
    show_default=True,
    help=u'The charge density cutoff energy in Ry.')

HUBBARD_U = OverridableOption(
    '-U',
    '--hubbard-u',
    nargs=2,
    multiple=True,
    type=click.Tuple([six.text_type, float]),
    help=u'Add a Hubbard U term to a specific kind.',
    metavar='<KIND MAGNITUDE>...')

HUBBARD_V = OverridableOption(
    '-V',
    '--hubbard-v',
    nargs=4,
    multiple=True,
    type=click.Tuple([int, int, int, float]),
    help=u'Add a Hubbard V interaction between two sites.',
    metavar='<SITE SITE TYPE MAGNITUDE>...')

HUBBARD_FILE = OverridableOption(
    '-H',
    '--hubbard-file',
    'hubbard_file_pk',
    type=types.DataParamType(sub_classes=('aiida.data:singlefile',)),
    help=u'SinglefileData containing Hubbard parameters from a HpCalculation to use as input for Hubbard V.')

STARTING_MAGNETIZATION = OverridableOption(
    '-M',
    '--starting-magnetization',
    nargs=2,
    multiple=True,
    type=click.Tuple([six.text_type, float]),
    help=u'Add a starting magnetization to a specific kind.',
    metavar='<KIND MAGNITUDE>...')

SMEARING = OverridableOption(
    '-S',
    '--smearing',
    nargs=2,
    default=(None, None),
    type=click.Tuple([six.text_type, float]),
    help=u'Add smeared occupations by specifying the type and amount of smearing.',
    metavar='<TYPE DEGAUSS>')
