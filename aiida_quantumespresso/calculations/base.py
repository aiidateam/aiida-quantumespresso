# -*- coding: utf-8 -*-
"""Defines a `CalcJob` base class for `aiida-quantumespresso`.

The custom `CalcJob` base class automatically sets the `invalidates_cache` attribute of exit codes based on the status
integer. All `CalcJob` implementations in `aiida-quantumespresso` must use this base class, not `aiida.engine.CalcJob`.
"""
from __future__ import absolute_import

from aiida.engine import CalcJob as _BaseCalcJob
from aiida.engine.processes.process_spec import CalcJobProcessSpec as _BaseCalcJobProcessSpec

__all__ = ('CalcJob',)


class CalcJobProcessSpec(_BaseCalcJobProcessSpec):
    """Process spec for `aiida-quantumespresso` `CalcJob` classes.

    Automatically sets the `invalidates_cache` flag to `True` if `exit_status < 400`, unless explicitly overridden.
    """

    def exit_code(self, status, label, message, invalidates_cache=None):
        """Add an exit code to the ProcessSpec.

        .. note: for any status smaller than `400` the `invalidates_cache` is automatically set to `True`, unless
            explicitly overridden.

        :param status: the exit status integer
        :param label: a label by which the exit code can be addressed
        :param message: a more detailed description of the exit code
        :param invalidates_cache: when set to `True`, a process exiting with this exit code will not be considered for
            caching
        """
        if invalidates_cache is None:
            invalidates_cache = (isinstance(status, int) and status < 400)

        try:
            super(CalcJobProcessSpec,
                  self).exit_code(status=status, label=label, message=message, invalidates_cache=invalidates_cache)
        except TypeError:  # For AiiDA version <1.1
            super(CalcJobProcessSpec, self).exit_code(status=status, message=message, label=label)


class CalcJob(_BaseCalcJob):  # pylint: disable=abstract-method
    """Custom `CalcJob` class for `aiida-quantumespresso` calculations."""

    _spec_class = CalcJobProcessSpec

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(CalcJob, cls).define(spec)
        spec.inputs['metadata']['options']['resources'].default = lambda: {'num_machines': 1}
