# -*- coding: utf-8 -*-
"""Workchain to perform a NEB calculation with climbing image."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, append_

from aiida_quantumespresso.common.types import ElectronicType, SpinType

from ...workflows.pw.base import PwBaseWorkChain
from .base import NebBaseWorkChain


class NebCiWorkChain(WorkChain):
    """Workchain to perform a NEB calculation with climbing image."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(NebBaseWorkChain, namespace='base',
            namespace_options={'help': 'Inputs for the `NebBaseWorkChain` for the main relax loop.'})

        spec.inputs.validator = cls.validate_inputs

        spec.outline(
            cls.run_noci_neb,
            cls.inspect_noci_neb,
            cls.run_ci_neb,
            cls.inspect_ci_neb,
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_NOCI',
            message='No-Climbing Image NebBaseWorkChain was excepted or killed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_CI',
            message='Climbing Image NebBaseWorkChain was excepted or killed')
        spec.expose_outputs(NebBaseWorkChain)
        # yapf: enable

    @classmethod
    def validate_inputs(cls, inputs, _):
        """Validate the top-level inputs."""
        if 'CLIMBING_IMAGES' in inputs['base']['neb']['parameters']:
            if inputs['base']['neb']['parameters']['PATH']['CI_scheme'] == 'auto':
                return 'Incopatible parameters: when using `CLIMBING_IMAGES`, `CI_scheme` cannot be `auto`.'
            inputs['base']['neb']['parameters']['PATH']['CI_scheme'] = 'manual'
        else:
            inputs['base']['neb']['parameters']['PATH']['CI_scheme'] = 'auto'

        cls.inputs = AttributeDict(inputs)

    @classmethod
    def get_builder_from_protocol(
        cls,
        code,
        images,
        protocol=None,
        overrides=None,
        electronic_type=ElectronicType.METAL,
        spin_type=SpinType.NONE,
        initial_magnetic_moments=None,
        options=None,
        **kargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param images: the ``TrajectoryData`` instance to use for initial guess of NEB images.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param electronic_type: indicate the electronic character of the system through ``ElectronicType`` instance.
        :param spin_type: indicate the spin polarization type to use through a ``SpinType`` instance.
        :param initial_magnetic_moments: optional dictionary that maps the initial magnetic moment of each kind to a
            desired value for a spin polarized calculation. Note that in case the ``starting_magnetization`` is also
            provided in the ``overrides``, this takes precedence over the values provided here. In case neither is
            provided and ``spin_type == SpinType.COLLINEAR``, an initial guess for the magnetic moments is used.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        pw_base = PwBaseWorkChain.get_builder_from_protocol(
            code,
            images.get_step_structure(-1),
            protocol=protocol,
            overrides=overrides,
            electronic_type=electronic_type,
            spin_type=spin_type,
            initial_magnetic_moments=initial_magnetic_moments,
            options=options,
            **kargs
        )
        #pylint: disable=no-member
        builder = cls.get_builder()
        builder.base.neb.code = code
        builder.base.neb.images = images
        builder.base.neb.pw.pseudos = pw_base.pw.pseudos
        builder.base.neb.pw.parameters = pw_base.pw.parameters
        builder.base.neb.metadata.options = pw_base.pw.metadata.options

        if 'kpoints' in pw_base:
            builder.base.kpoints = pw_base['kpoints']
        else:
            builder.base.kpoints_distance = orm.Float(pw_base['kpoints_distance'])
        builder.base.kpoints_force_parity = orm.Bool(pw_base['kpoints_force_parity'])
        builder.base.max_iterations = orm.Int(pw_base['max_iterations'])
        # pylint: enable=no-member
        return builder

    def run_noci_neb(self):
        """Run the `NebBaseWorkChain` to run a `NebCalculation`."""

        inputs = AttributeDict(self.exposed_inputs(NebBaseWorkChain, 'base'))

        parameters = inputs.neb.parameters.get_dict()
        parameters['PATH']['CI_scheme'] = 'no-CI'
        inputs.neb.parameters = parameters

        # Set the `CALL` link label
        inputs.metadata.call_link_label = 'no_climbing_image'

        running = self.submit(NebBaseWorkChain, **inputs)
        self.report(f'launching NebBaseWorkChain<{running.pk}>')

        return ToContext(workchains=append_(running))

    def inspect_noci_neb(self):
        """Inspect the results of the no-CI `NebBaseWorkChain`."""
        workchain = self.ctx.workchains[-1]

        if workchain.exit_status != 0:
            self.report('No-Climbing Image NebBaseWorkChain was excepted or killed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_NOCI

        self.ctx.remote_folder = workchain.outputs.remote_folder
        return

    def run_ci_neb(self):
        """Run a Climbing Image `NebBaseWorkChain`."""

        inputs = AttributeDict(self.exposed_inputs(NebBaseWorkChain, 'base'))
        inputs.neb.parent_folder = self.ctx.remote_folder
        parameters = inputs.neb.parameters.get_dict()
        parameters['PATH']['restart_mode'] = 'restart'
        inputs.neb.parameters = parameters
        # Set the `CALL` link label
        inputs.metadata.call_link_label = 'climbing_image'
        running = self.submit(NebBaseWorkChain, **inputs)

        self.report(f'launching NebBaseWorkChain<{running.pk}>')

        return ToContext(workchains=append_(running))

    def inspect_ci_neb(self):
        """Inspect the results of the Climbing Image `NebBaseWorkChain`."""
        workchain = self.ctx.workchains[-1]

        if workchain.exit_status != 0:
            self.report('Climbing Image NebBaseWorkChain was excepted or killed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_CI
        self.out_many(self.exposed_outputs(workchain, NebBaseWorkChain))

        return
