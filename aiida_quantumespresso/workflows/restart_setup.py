# -*- coding: utf-8 -*-
##############################################################################################################
"""Workchain to generate a restart remote folder for a Quantum ESPRESSO calculation."""

from aiida.orm import RemoteData, FolderData
from aiida.engine import WorkChain, ToContext
from aiida.plugins import WorkflowFactory, CalculationFactory

from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin
from aiida_quantumespresso.utils.transfer import get_transfer_builder

TransferCalcjob = CalculationFactory('core.transfer')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')

# pylint: disable=f-string-without-interpolation


##############################################################################################################
class RestartSetupWorkChain(ProtocolMixin, WorkChain):
    """Workchain to generate a restart remote folder for a Quantum ESPRESSO calculation.

    It consists of two steps:

        1. TransferCalcjob: takes the content of a ``FolderData`` node and copies it into the remote computer
           into a ``RemoteData`` folder with the correct folder structure. The original ``FolderData`` needs
           to have the necessary files in the right internal path (check the ``get_transfer_builder`` utility
           function and/or use it to retrieve the densities to have this already taken care of)

        2. PwBaseWorkChain: it runs an NSCF calculation using the previously created ``RemoteData`` as its
           ``parent_folder``. This will re-generate all the wavefunctions in the running directory, which
           are necessary for launching any other kind of QE calculation in restart mode.

    """

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)

        spec.expose_inputs(
            TransferCalcjob,
            namespace='transfer',
            namespace_options={
                'validator': validate_transfer,
                'populate_defaults': False,
                'help': 'Inputs for the `TransferCalcjob` to put the data on the cluster.',
            }
        )

        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='nscf',
            exclude=('clean_workdir', 'pw.parent_folder'),
            namespace_options={
                'validator': validate_nscf,
                'populate_defaults': False,
                'help': 'Inputs for the `PwBaseWorkChain` for the NSCF calculation.',
            }
        )

        spec.inputs.validator = validate_inputs

        spec.outline(
            cls.run_transfer,
            cls.inspect_transfer,
            cls.run_nscf,
            cls.inspect_nscf,
            cls.results,
        )

        spec.output(
            'remote_data',
            valid_type=RemoteData,
            help='The output node with the folder to be used as parent folder for other calculations.',
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_TRANSFER', message='The TransferCalcjob sub process failed.')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_NSCF', message='The ncf PwBasexWorkChain sub process failed.')

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        raise NotImplementedError(f'`get_protocol_filepath` method not yet implemented in `RestartSetupWorkChain`')

    @classmethod
    def get_builder_from_protocol(cls, folder_data, structure, code, protocol=None, overrides=None, **kwargs):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param data_source: the ``FolderData`` node containing the density (and the rest of the restart data).
        :param structure: the ``StructureData`` instance required to run the NSF calculation.
        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin, required to
            run the NSF calculation.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol``
            of all the sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        # inputs = cls.get_protocol_inputs(protocol, overrides)

        builder = cls.get_builder()

        track = kwargs.get('track', False)
        transfer = get_transfer_builder(folder_data, computer=code.computer, track=track)
        # The cls.builder seems to be checking that the metadata.options of each calcjob has
        # a resources entry? There should be an exception to this for the CalcJobs that skip
        # the submission, but until that is implemented in a supported version of aiida-core
        # this patch needs to remain here:
        transfer['metadata']['options']['resources'] = {}
        builder.transfer = transfer

        nscf_args = (code, structure, protocol)
        nscf_kwargs = kwargs
        nscf_kwargs['overrides'] = {}
        if overrides is not None:
            nscf_kwargs['overrides'] = overrides.get('nscf', None)

        # This is for easily setting defaults at each level of:
        # [overrides.nscf].pw.parameters.CONTROL.calculation
        last_layer = nscf_kwargs['overrides']
        last_layer = last_layer.setdefault('pw', {})
        last_layer = last_layer.setdefault('parameters', {})
        last_layer = last_layer.setdefault('CONTROL', {})

        if last_layer.setdefault('calculation', 'nscf') != 'nscf':
            bad_value = last_layer['calculation']
            raise ValueError(
                f'The internal PwBaseWorkChain is for running an NSCF calculation, '
                f'this should not be overriden. '
                f'(Found overrides.nscf.pw.parameters.CONTROL.calculation=`{bad_value}`)'
            )

        nscf = PwBaseWorkChain.get_builder_from_protocol(*nscf_args, **nscf_kwargs)
        nscf['pw'].pop('parent_folder', None)
        nscf.pop('clean_workdir', None)
        builder.nscf = nscf

        return builder

    def run_transfer(self):
        """Run the TransferCalcjob to put the data in the remote computer."""
        inputs = self.exposed_inputs(TransferCalcjob, namespace='transfer')
        running = self.submit(TransferCalcjob, **inputs)
        self.report(f'launching TransferCalcjob<{running.pk}> for put the data into the remote computer')
        return ToContext(transfer_calcjob=running)

    def inspect_transfer(self):
        """Verify that the TransferCalcjob to get data finished successfully."""
        calcjob_node = self.ctx.transfer_calcjob

        if not calcjob_node.is_finished_ok:
            self.report(f'TransferCalcjob failed with exit status {calcjob_node.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_TRANSFER

        self.ctx.remote_parent = calcjob_node.outputs.remote_folder

    def run_nscf(self):
        """Run the PwBaseWorkChain in nscf mode on the restart folder."""
        inputs = self.exposed_inputs(PwBaseWorkChain, namespace='nscf')
        inputs['metadata']['call_link_label'] = 'nscf'
        inputs['pw']['parent_folder'] = self.ctx.remote_parent

        running = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{running.pk}> in nscf mode')
        return ToContext(workchain_nscf=running)

    def inspect_nscf(self):
        """Verify that the PwBaseWorkChain for the scf run finished successfully."""
        workchain = self.ctx.workchain_nscf

        if not workchain.is_finished_ok:
            self.report(f'scf PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_NSCF

        self.ctx.remote_data = workchain.outputs.remote_folder

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain succesfully completed')
        self.out('remote_data', self.ctx.remote_data)


##############################################################################################################
def validate_transfer(value, _):
    """Validate the inputs of the transfer input namespace."""

    # Check that the source node is there and is of right type
    if 'source_nodes' not in value:
        return f'The inputs of the transfer namespace were not set correctly: {value}'

    source_nodes = value['source_nodes']
    if 'source_node' not in source_nodes:
        return f'The `source_nodes` in the transfer namespace was not set correctly: {source_nodes}'

    source_node = source_nodes['source_node']
    if not isinstance(source_node, FolderData):
        return f'The `source_node` in the transfer namespace is not `FolderData`: {source_node}'

    # Check that the files are in the source node
    error_message = ''
    if 'data-file-schema.xml' not in source_node.list_object_names():
        error_message += f'Missing `data-file-schema.xml` on node PK={source_node.pk}\n'

    if 'charge-density.dat' not in source_node.list_object_names():
        error_message += f'Missing `charge-density.dat` on node PK={source_node.pk}\n'

    if len(error_message) > 0:
        return error_message


def validate_nscf(value, _):
    """Validate the inputs of the nscf input namespace."""
    parameters = value['pw']['parameters'].get_dict()
    if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'nscf':
        return '`CONTOL.calculation` in `nscf.pw.parameters` is not set to `nscf`.'


def validate_inputs(inputs, _):
    """Validate the inputs of the entire input namespace."""
    computer_transfer = inputs['transfer']['metadata']['computer']
    computer_nscf = inputs['nscf']['pw']['code'].computer

    if computer_transfer.pk != computer_nscf.pk:
        return (
            f'The computer where the files are being copied ({computer_transfer}) '
            f'is not where the code resides ({computer_nscf})'
        )
