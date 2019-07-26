# -*- coding: utf-8 -*-
"""Command line scripts to launch a `NebCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..cli import calculation_launch
from ..utils import launch
from ..utils import options as options_qe
from ..utils import validate

@calculation_launch.command('neb')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.neb'))
@click.option(
    '-s',
    '--structures',
    nargs=2,
    type=types.DataParamType(sub_classes=('aiida.data:structure',)),
    help='Two StructureData nodes, the initial and final structures',
    required=True)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.STARTING_MAGNETIZATION()
@options_qe.SMEARING()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, structures, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, starting_magnetization,
                       smearing, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """
    Run a NebCalculation.
    Note that some parameters are hardcoded.
    """
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    # TODO: move some hardcoded parameters to unit tests
    # TODO: because I cannot give a structure-specific settings dict, I can't apply these constraints only to one of the
    #  two boundary images. Do I need to?
    settings = {
        'fixed_coords':
            [[False, True, True],
             [True, True, True],
             [False, True, True]],
    }

    pw_parameters = {
        'CONTROL': {
            'calculation': 'relax'
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
            'occupations': 'smearing',
            'degauss': 0.003,
            'nspin': 2,
            'starting_magnetization': 0.5,
        },
        'ELECTRONS': {
            'conv_thr': 1e-8,
            'mixing_beta': 0.3,
        }
    }

    neb_parameters = {
        'PATH': {
            'nstep_path': 20,
            'ds': 2.,
            'opt_scheme': 'broyden',
            'num_of_images': 7,
            'k_max': 0.3,
            'k_min': 0.2,
            'CI_scheme': "auto",
            'path_thr': 0.1,
        }
    }

    for structure in structures:
        try:
            validate.validate_starting_magnetization(structure, pw_parameters, starting_magnetization)
        except ValueError as exception:
            raise click.BadParameter(str(exception))

        try:
            validate.validate_smearing(pw_parameters, smearing)
        except ValueError as exception:
            raise click.BadParameter(str(exception))

    inputs = {
        'code': code,
        'first_structure': structures[0],
        'last_structure': structures[1],
        'pw': {
            'pseudos': get_pseudos_from_structure(structures[0], pseudo_family),
            'kpoints': kpoints_mesh,
            'parameters': Dict(dict=pw_parameters),
        },
        'parameters': Dict(dict=neb_parameters),
        'settings': Dict(dict=settings),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.neb'), daemon, **inputs)
