# -*- coding: utf-8 -*-
"""TCOD export plugin for `CpCalculations`."""
try:
    from aiida_tcod.tools.dbexporters.tcod import BaseTcodtranslator  # pylint: disable=import-error
except ImportError as exception:
    raise ImportError('dependency `aiida-tcod` not installed; run `pip install aiida-tcod` to do so.') from exception


class CpTcodtranslator(BaseTcodtranslator):
    """TCOD export plugin for `CpCalculations`."""

    # pylint: disable=abstract-method

    _plugin_type_string = 'quantumespresso.cp.CpCalculation'

    @classmethod
    def get_software_package(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the package or program name that was used to produce the structure.

        Only package or program name should be used, e.g. 'VASP', 'psi3', 'Abinit', etc.
        """
        return 'Quantum ESPRESSO'

    @classmethod
    def get_number_of_electrons(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the number of electrons."""
        parameters = calc.outputs.output_parameters
        if 'number_of_electrons' not in parameters.attrs():
            return None
        return parameters.get_attr('number_of_electrons')

    @classmethod
    def get_computation_wallclock_time(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the computation wallclock time in seconds."""
        parameters = calc.outputs.output_parameters
        if 'wall_time_seconds' not in parameters.attrs():
            return None
        return parameters.get_attr('wall_time_seconds')
