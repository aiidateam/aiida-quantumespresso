# -*- coding: utf-8 -*-
##############################################################################################################
"""Return process builders ready for transferring Quantum ESPRESSO density restart files."""

import warnings

from aiida.orm import FolderData, RemoteData, Dict
from aiida.engine import calcfunction
from aiida.plugins import CalculationFactory


def get_transfer_builder(data_source, computer=None, track=False):
    """Create a `ProcessBuilder` for `TransferCalcjob`from a data_source.

    The data_source can be of either `RemoteData` or `FolderData`:

      - `RemoteData`: generate a set of instructions so that the density restart data will be taken from the
        remote computer specified by the node and into the local aiida DB.

      - `FolderData`: generate a set of instructions so that the density restart data will be taken from the
        local aiida DB and into the provided computer (which has to be given as an extra parameter).

    :param data_source: the node instance from which to take the density data
    :param computer: if `data_source` is a `FolderData` node, the remote computer to which to transfer the data must
        be specified here
    :param track: boolean, if True, the generation of the instructions will be done through a calcfunction from the
        data_source as input (and thus be tracked as such in the provenance)
    :return: a `ProcessBuilder` instance configured for launching a `TransferCalcjob`
    """
    builder = CalculationFactory('core.transfer').get_builder()
    builder.source_nodes = {'source_node': data_source}

    if isinstance(data_source, FolderData):
        if computer is None:
            raise ValueError('No computer was provided for setting up a transfer to a remote.')
        builder.metadata['computer'] = computer

    elif isinstance(data_source, RemoteData):
        if computer is not None:
            warnings.warn(
                f'Computer `{computer}` provided will be ignored '
                f'(using `{data_source.computer}` from the RemoteData input `{data_source}`)'
            )
        builder.metadata['computer'] = data_source.computer

    if track:
        builder.instructions = generate_instructions(data_source)['instructions']
    else:
        builder.instructions = generate_instructions_untracked(data_source)

    return builder


##############################################################################################################
def generate_instructions_untracked(source_folder):
    """Generate the instruction node to be used for copying the files."""

    # Paths in the QE run folder
    schema_qepath = 'out/aiida.save/data-file-schema.xml'
    charge_qepath = 'out/aiida.save/charge-density.dat'
    pawtxt_qepath = 'out/aiida.save/paw.txt'

    # Paths in the local node
    schema_dbpath = 'data-file-schema.xml'
    charge_dbpath = 'charge-density.dat'
    pawtxt_dbpath = 'paw.txt'

    # Transfer from local to remote
    if isinstance(source_folder, FolderData):
        instructions = {'retrieve_files': False, 'local_files': []}
        instructions['local_files'].append(('source_node', schema_dbpath, schema_qepath))
        instructions['local_files'].append(('source_node', charge_dbpath, charge_qepath))

        if 'paw.txt' in source_folder.list_object_names():
            instructions['local_files'].append(('source_node', pawtxt_dbpath, pawtxt_qepath))

    # Transfer from remote to local
    elif isinstance(source_folder, RemoteData):
        instructions = {'retrieve_files': True, 'symlink_files': []}
        instructions['symlink_files'].append(('source_node', schema_qepath, schema_dbpath))
        instructions['symlink_files'].append(('source_node', charge_qepath, charge_dbpath))
        instructions['symlink_files'].append(('source_node', pawtxt_qepath, pawtxt_dbpath))

    return Dict(dict=instructions)


@calcfunction
def generate_instructions(source_folder):
    """Auxiliary function to keep provenance track of the generation of the instructions."""
    output_node = generate_instructions_untracked(source_folder)
    return {'instructions': output_node}
