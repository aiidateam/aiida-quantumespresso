# -*- coding: utf-8 -*-
"""Command line scripts to launch a `EpwCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..utils import launch
from ..utils import options as options_qe
from . import cmd_launch


@cmd_launch.command('epw')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.epw'))
@click.option('--pw-parent', type=types.CalculationParamType(), help='The parent pw.x calculation.')
@click.option('--ph-parent', type=types.CalculationParamType(), help='The parent pw.x calculation.')
@options_qe.KPOINTS_MESH(default=[6, 6, 6])
@options_qe.QPOINTS_MESH(default=[2, 2, 2])
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, kpoints_mesh, qpoints_mesh, pw_parent, ph_parent, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a EpwCalculation."""
    from aiida import orm
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    # Check that the parent calculation node comes from quantumespresso.pw and quantumespresso.ph
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.pw'
    if pw_parent.process_type != expected_process_type:
        raise click.BadParameter(
            'The input calculation node has a process_type: {}; should be {}'.format(
                pw_parent.process_type, expected_process_type
            )
        )

    pw_parent_folder = pw_parent.get_outgoing(node_class=orm.RemoteData, link_label_filter='remote_folder').one().node

    expected_process_type = 'aiida.calculations:quantumespresso.ph'
    if ph_parent.process_type != expected_process_type:
        raise click.BadParameter(
            'The input calculation node has a process_type: {}; should be {}'.format(
                ph_parent.process_type, expected_process_type
            )
        )

    ph_parent_folder = ph_parent.get_outgoing(node_class=orm.RemoteData, link_label_filter='remote_folder').one().node


    inputs = {
        'code': code,
        'qpoints': qpoints_mesh,
        'kpoints': kpoints_mesh,
        'parameters': orm.Dict(dict={'INPUTEPW': {}}),
        'parent_folder': pw_parent_folder,
        'parent_folder_ph': ph_parent_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.epw'), daemon, **inputs)
