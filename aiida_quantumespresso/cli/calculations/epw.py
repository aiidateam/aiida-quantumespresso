# -*- coding: utf-8 -*-
"""Command line scripts to launch a `EpwCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click
import numpy as np

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..utils import launch
from ..utils import options as options_qe
from . import cmd_launch


@cmd_launch.command('epw')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.epw'))
@click.option('--pw-nscf-parent', type=types.CalculationParamType(), help='The parent pw.x nscf calculation.')
@click.option('--ph-parent', type=types.CalculationParamType(), help='The parent ph.x calculation.')
@options_qe.KPOINTS_MESH(default=[6, 6, 6])
@options_qe.QPOINTS_MESH(default=[2, 2, 2])
@options_qe.QIBZ(default=((0.0, 0.0, 0.0),
    (0.353553391, -0.353553391, -0.353553391), (0.0, 0.0, -0.707106781)))
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()

def launch_calculation(code, kpoints_mesh, qpoints_mesh, qibz, pw_nscf_parent,
        ph_parent, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a EpwCalculation."""
    from aiida import orm
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    qibz_node = orm.ArrayData()
    qibz_node.set_array('qibz', np.array(qibz))
    # Check that the parent calculation node comes from quantumespresso.pw and quantumespresso.ph
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.pw'
    if pw_nscf_parent.process_type != expected_process_type:
        raise click.BadParameter(
            'The input calculation node has a process_type: {}; should be {}'.format(
                pw_nscf_parent.process_type, expected_process_type
            )
        )

    pw_nscf_parent_folder = pw_nscf_parent.get_outgoing(node_class=orm.RemoteData,
            link_label_filter='remote_folder').one().node

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
        'qibz': qibz_node,
        'parameters': orm.Dict(dict={'INPUTEPW': {
                                         'nbndsub'     : 8,
                                         'elph'        : True,
                                         'kmaps'       : False,
                                         'epbwrite'    : True,
                                         'epbread'     : False,
                                         'epwwrite'    : True,
                                         'epwread'     : False,
                                         'proj(1)'     : 'Si : sp3',
                                         'elecselfen'  : True,
                                         'nkf1'        : 6,
                                         'nkf2'        : 6,
                                         'nkf3'        : 6,
                                         'nqf1'        : 2,
                                         'nqf2'        : 2,
                                         'nqf3'        : 2,
              				 'wannierize'  : True,
                                         'dvscf_dir'   : './save/',
                                         'dis_win_max' : 18,
                                         'dis_froz_max': 8.5
                                         }}),
        'parent_folder_nscf': pw_nscf_parent_folder,
        'parent_folder_ph': ph_parent_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.epw'), daemon, **inputs)

# Example to run an EPW calculation with AiiDA
#
# 1) Configure the pw.x, ph.x and epw.x codes called here pw@localhost, ph@localhost and epw@localhost
#
# 2) Run scf calculation
# aiida-quantumespresso calculation launch pw -X pw@localhost -p   SSSP_1.1_efficiency
# Retrieve the node PK (named here NODE_PK_SCF)
#
# 3) Run a phonon calculation
# aiida-quantumespresso calculation launch ph -X ph@localhost  -C NODE_PK_SCF
# Retrieve the node PK (named here NODE_PK_PH)
#
# 4) Run an nscf calculation
# This is not standard and needs to be done within the verdi shell:
# import os
# import numpy as np
# from aiida.engine import submit
# from aiida import orm
#
# PwCalculation = CalculationFactory('quantumespresso.pw')
#
# first_pw = load_node(2540)
# builder = first_pw.get_builder_restart()
# updated_parameters = builder.parameters.get_dict()
# updated_parameters['CONTROL']['calculation'] = 'nscf'
# updated_parameters['SYSTEM']['nbnd'] = 10
#
# KpointsData = DataFactory('array.kpoints')
# kpoints = KpointsData()
#
# klist = np.zeros((216, 3))
# tt = 0
# for ii in np.arange(0, 1, 1.0/6):
#   for jj in np.arange(0, 1, 1.0/6):
#     for kk in np.arange(0, 1, 1.0/6):
#       klist[tt, :] = np.array([ii, jj, kk])
#       tt += 1
# kpoints.set_kpoints(klist, cartesian = False, weights= np.ones(216)*1.0/216)
# kpoints.store()
#
# builder.kpoints = kpoints
# builder.parameters = Dict(dict=updated_parameters)
#
# builder.parent_folder = first_pw.outputs.remote_folder
#
# submit(builder)
#
# Record the PK number from that calculation (named here NODE_PK_NSCF)
#
# 5) Run an EPW calculation
# aiida-quantumespresso calculation launch epw -X epw@localhost --pw-nscf-parent NODE_PK_NSCF  --ph-parent NODE_PK_PH
#
# 6) Retrive your data from the EPW calculation
# verdi process list -a to get the NODE_PK_EPW from the EPW calculation
# verdi calcjob gotocomputer NODE_PK_EPW
