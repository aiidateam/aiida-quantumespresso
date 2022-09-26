# -*- coding: utf-8 -*-
"""Command line scripts to launch a `EpwCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators
import click

from . import cmd_launch
from ..utils import launch, options


@cmd_launch.command('epw')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.epw'))
@click.option('--pw-nscf-parent', type=types.CalculationParamType(), help='The parent pw.x nscf calculation.')
@click.option('--ph-parent', type=types.CalculationParamType(), help='The parent ph.x calculation.')
@options.KPOINTS_MESH(default=[6, 6, 6])
@options.QPOINTS_MESH(default=[2, 2, 2])
@options.KFPOINTS_MESH(default=[6, 6, 6])
@options.QFPOINTS_MESH(default=[2, 2, 2])
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_calculation(
    code, kpoints_mesh, qpoints_mesh, kfpoints_mesh, qfpoints_mesh, pw_nscf_parent, ph_parent, max_num_machines,
    max_wallclock_seconds, with_mpi, daemon
):
    """Run a EpwCalculation."""
    from aiida import orm
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    # Check that the parent calculation node comes from quantumespresso.pw and quantumespresso.ph
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.pw'
    if pw_nscf_parent.process_type != expected_process_type:
        raise click.BadParameter(
            f'--pw-nscf-parent node has process_type: {pw_nscf_parent.process_type}; should be {expected_process_type}'
        )

    pw_nscf_parent_folder = pw_nscf_parent.base.links.get_outgoing(
        node_class=orm.RemoteData, link_label_filter='remote_folder'
    ).one().node

    expected_process_type = 'aiida.calculations:quantumespresso.ph'
    if ph_parent.process_type != expected_process_type:
        raise click.BadParameter(
            f'--ph-parent has process_type: {ph_parent.process_type}; should be {expected_process_type}'
        )

    ph_parent_folder = ph_parent.base.links.get_outgoing(node_class=orm.RemoteData,
                                                         link_label_filter='remote_folder').one().node

    parameters = {
        'INPUTEPW': {
            'nbndsub': 8,
            'elph': True,  # default is false
            'epbwrite': True,  # default is false
            'epwwrite': True,  # default is false
            'proj(1)': 'Si : sp3',
            'elecselfen': True,
            'wannierize': True,
            'dvscf_dir': './save/',
            'dis_win_max': 18,
            'dis_froz_max': 8.5
        }
    }

    inputs = {
        'code': code,
        'qpoints': qpoints_mesh,
        'kpoints': kpoints_mesh,
        'qfpoints': qfpoints_mesh,
        'kfpoints': kfpoints_mesh,
        'parameters': orm.Dict(parameters),
        'parent_folder_nscf': pw_nscf_parent_folder,
        'parent_folder_ph': ph_parent_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.epw'), daemon, **inputs)
