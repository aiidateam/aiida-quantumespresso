# -*- coding: utf-8 -*-
"""
Immigrate a Quantum Espresso pw.x job that was not run using AiiDa.
"""

from __future__ import absolute_import

from aiida.common.datastructures import CalcJobState
from aiida.common.folders import SandboxFolder
from aiida.common.links import LinkType
from aiida.engine import ProcessState
from aiida.engine.daemon.execmanager import retrieve_calculation
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager


def immigrate_existing(builder, remote_data):
    """Immigrate a Calculation that was not run using AiiDa.

    :param builder: a populated builder instance for a CalcJob
    :type builder: aiida.engine.processes.builder.ProcessBuilder
    :param remote_data: a remote data folder, containing the output files required for parsing
    :type remote_data: aiida.orm.RemoteData

    :rtype: aiida.orm.CalcJobNode

    """
    # initialise calcjob
    runner = get_manager().get_runner()
    pw_calc_cls = builder._process_class
    process = instantiate_process(runner, pw_calc_cls, **builder)
    calc_node = process.node

    # prepare for submission
    with SandboxFolder() as temp_folder:
        calc_info = process.presubmit(temp_folder)
        calc_node.put_object_from_tree(temp_folder.abspath, force=True)

    # link remote folder to calc_node
    if not remote_data.is_stored:
        remote_data.store()
    remote_data.add_incoming(calc_node, link_type=LinkType.CREATE, link_label='remote_folder')
    calc_node.set_remote_workdir(remote_data.get_remote_path())
    transport = remote_data.computer.get_transport()

    with SandboxFolder() as temp_retrieved:
        # retrieved output files
        retrieve_calculation(calc_node, transport, temp_retrieved.abspath)
        # parse output
        calc_node.set_state(CalcJobState.PARSING)
        exit_code = process.parse(temp_retrieved.abspath)
    # link outgoing nodes
    process.update_outputs()

    # finalise calc node
    calc_node.delete_state()
    calc_node.delete_checkpoint()
    calc_node.set_process_state(ProcessState.FINISHED)
    calc_node.set_exit_status(exit_code.status)
    calc_node.set_exit_message(exit_code.message)
    calc_node.seal()

    # record that the node was created via immigration
    calc_node.set_extra('immigrated', True)
    calc_node.set_extra('immigration_func', __name__)

    return calc_node
