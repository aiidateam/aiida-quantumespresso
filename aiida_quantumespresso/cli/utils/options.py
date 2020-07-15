# -*- coding: utf-8 -*-
"""Pre-defined overridable options for commonly used command line interface parameters."""
import click

from aiida.cmdline.params import types
from aiida.cmdline.params.options import OverridableOption

from . import validate

STRUCTURE = OverridableOption(
    '-s', '--structure', type=types.DataParamType(sub_classes=('aiida.data:structure',)), help='StructureData node.'
)

PSEUDO_FAMILY = OverridableOption(
    '-p', '--pseudo-family', 'pseudo_family', type=click.STRING, help='Pseudo potential family name.'
)

KPOINTS_DISTANCE = OverridableOption(
    '-K',
    '--kpoints-distance',
    type=click.FLOAT,
    default=0.5,
    show_default=True,
    help='The minimal distance between k-points in reciprocal space in inverse Ångström.'
)

KPOINTS_MESH = OverridableOption(
    '-k',
    '--kpoints-mesh',
    'kpoints_mesh',
    nargs=3,
    type=click.INT,
    show_default=True,
    callback=validate.validate_kpoints_mesh,
    help='The number of points in the kpoint mesh along each basis vector.'
)

QPOINTS_MESH = OverridableOption(
    '-q',
    '--qpoints-mesh',
    'qpoints_mesh',
    nargs=3,
    type=click.INT,
    show_default=True,
    callback=validate.validate_kpoints_mesh,
    help='The number of points in the qpoint mesh along each basis vector.'
)

KFPOINTS_MESH = OverridableOption(
    '-kf',
    '--kfpoints-mesh',
    'kfpoints_mesh',
    nargs=3,
    type=click.INT,
    show_default=True,
    callback=validate.validate_kpoints_mesh,
    help='The number of points in the fine kpoint mesh along each basis vector.'
)

QFPOINTS_MESH = OverridableOption(
    '-qf',
    '--qfpoints-mesh',
    'qfpoints_mesh',
    nargs=3,
    type=click.INT,
    show_default=True,
    callback=validate.validate_kpoints_mesh,
    help='The number of points in the fine qpoint mesh along each basis vector.'
)

QIBZ = OverridableOption(
    '--qpoint-ibz',
    'qibz',
    nargs=3,
    multiple=True,
    type=click.FLOAT,
    show_default=True,
    help='The IBZ q-point list. Must be the same as the previous PH calculation.'
)

MAX_NUM_MACHINES = OverridableOption(
    '-m',
    '--max-num-machines',
    type=click.INT,
    default=1,
    show_default=True,
    help='The maximum number of machines (nodes) to use for the calculations.'
)

MAX_WALLCLOCK_SECONDS = OverridableOption(
    '-w',
    '--max-wallclock-seconds',
    type=click.INT,
    default=1800,
    show_default=True,
    help='the maximum wallclock time in seconds to set for the calculations.'
)

WITH_MPI = OverridableOption(
    '-i', '--with-mpi', is_flag=True, default=False, show_default=True, help='Run the calculations with MPI enabled.'
)

PARENT_FOLDER = OverridableOption(
    '-P',
    '--parent-folder',
    'parent_folder',
    type=types.DataParamType(sub_classes=('aiida.data:remote',)),
    show_default=True,
    required=False,
    help='The PK of a parent remote folder (for restarts).'
)

DAEMON = OverridableOption(
    '-d',
    '--daemon',
    is_flag=True,
    default=False,
    show_default=True,
    help='Submit the process to the daemon instead of running it locally.'
)

AUTOMATIC_PARALLELIZATION = OverridableOption(
    '-a',
    '--automatic-parallelization',
    is_flag=True,
    default=False,
    show_default=True,
    help='Enable the automatic parallelization option of the workchain.'
)

CLEAN_WORKDIR = OverridableOption(
    '-x',
    '--clean-workdir',
    is_flag=True,
    default=False,
    show_default=True,
    help='Clean the remote folder of all the launched calculations after completion of the workchain.'
)

ECUTWFC = OverridableOption(
    '-W', '--ecutwfc', type=click.FLOAT, default=30., show_default=True, help='The plane wave cutoff energy in Ry.'
)

ECUTRHO = OverridableOption(
    '-R',
    '--ecutrho',
    type=click.FLOAT,
    default=240.,
    show_default=True,
    help='The charge density cutoff energy in Ry.'
)

HUBBARD_U = OverridableOption(
    '-U',
    '--hubbard-u',
    nargs=2,
    multiple=True,
    type=click.Tuple([str, float]),
    help='Add a Hubbard U term to a specific kind.',
    metavar='<KIND MAGNITUDE>...'
)

HUBBARD_V = OverridableOption(
    '-V',
    '--hubbard-v',
    nargs=4,
    multiple=True,
    type=click.Tuple([int, int, int, float]),
    help='Add a Hubbard V interaction between two sites.',
    metavar='<SITE SITE TYPE MAGNITUDE>...'
)

HUBBARD_FILE = OverridableOption(
    '-H',
    '--hubbard-file',
    'hubbard_file_pk',
    type=types.DataParamType(sub_classes=('aiida.data:singlefile',)),
    help='SinglefileData containing Hubbard parameters from a HpCalculation to use as input for Hubbard V.'
)

STARTING_MAGNETIZATION = OverridableOption(
    '-M',
    '--starting-magnetization',
    nargs=2,
    multiple=True,
    type=click.Tuple([str, float]),
    help='Add a starting magnetization to a specific kind.',
    metavar='<KIND MAGNITUDE>...'
)

SMEARING = OverridableOption(
    '-S',
    '--smearing',
    nargs=2,
    default=(None, None),
    type=click.Tuple([str, float]),
    help='Add smeared occupations by specifying the type and amount of smearing.',
    metavar='<TYPE DEGAUSS>'
)
