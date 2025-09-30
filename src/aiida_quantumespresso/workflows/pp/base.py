"""Workchain to run a Quantum ESPRESSO pp.x calculation with automated error handling and restarts."""

from __future__ import annotations

from aiida import orm
from aiida.common import AttributeDict
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, while_

from aiida_quantumespresso.calculations.pp import PpCalculation
from aiida_quantumespresso.common.types import PostProcessQuantity


class PpBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO pp.x calculation with automated error handling and restarts."""

    _process_class = PpCalculation

    PP_CONFIG: dict[str, dict] = {
        'charge_density': {
            'plot_num': 0,
        },
        'potential': {
            'plot_num': 11,
        },
    }

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        spec.expose_inputs(PpCalculation, namespace='pp')
        spec.outline(
            cls.setup,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.expose_outputs(PpCalculation, exclude=('output_data_multiple',))

        spec.exit_code(
            206,
            'ERROR_SUB_PROCESS_FAILED',
            message='The PP calculation failed.',
        )

    @classmethod
    def get_parameters(cls, quantity: str) -> orm.Dict:
        """Return the parameters based on the requested quantity.

        :param quantity: The physical quantity to compute, as a string.
        :return: A `Dict` instance with the parameters.
        """
        config = cls.PP_CONFIG.get(quantity, {})
        parameters = {
            'INPUTPP': {
                'plot_num': config.get('plot_num'),
            },
            'PLOT': {
                'iflag': 3,
            },
        }
        return orm.Dict(parameters)

    @classmethod
    def get_builder_from_quantity(
        cls,
        code: orm.Code | str,
        quantity: PostProcessQuantity,
        parent_folder: orm.RemoteData,
        options: dict | None = None,
        **_,
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: The ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param quantity: The physical quantity to compute, as a ``PostProcessQuantity`` instance.
        :param parent_folder: The `RemoteData` node that contains the output of a completed `PwCalculation`.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :return: A process builder instance with all inputs defined ready for launch.
        :raises NotImplementedError: if the specified quantity is not supported.
        """
        if isinstance(code, str):
            code = orm.load_code(code)

        type_check(code, orm.AbstractCode)
        type_check(quantity, PostProcessQuantity)

        if quantity not in tuple(PostProcessQuantity):
            raise NotImplementedError(f'quantity `{quantity}` is not supported.')

        builder = cls.get_builder()
        builder.pp = {
            'code': code,
            'parent_folder': parent_folder,
            'parameters': cls.get_parameters(quantity.value),
            'metadata': {'options': options or {}},
        }

        return builder

    def setup(self) -> None:
        super().setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PpCalculation, namespace='pp'))
