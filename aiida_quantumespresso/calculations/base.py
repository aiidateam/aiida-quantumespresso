# -*- coding: utf-8 -*-
"""Defines a CalcJob base class for aiida-quantumespresso.

The custom CalcJob base class automatically sets the `invalidates_cache`
attribute of exit codes based on the status integer.
"""

from __future__ import absolute_import

from aiida.engine import ExitCode, CalcJob as _BaseCalcJob
from aiida.engine.processes.process_spec import CalcJobProcessSpec as _BaseCalcJobProcessSpec

__all__ = ('CalcJob',)


class CalcJobProcessSpec(_BaseCalcJobProcessSpec):
    """Process spec for aiida-quantumespresso CalcJob classes.

    Automatically sets the `invalidates_cache` flag to `True` for an
    exit status smaller than 400, unless explicitly overriden.
    """

    def exit_code(self, status, label, message, invalidates_cache=None):
        """Add an exit code to the ProcessSpec.

        .. note: for any status smaller than `400` the `invalidates_cache` is
            automatically set to `True`, unless explicitly overridden.

        :param status: the exit status integer
        :param label: a label by which the exit code can be addressed
        :param message: a more detailed description of the exit code
        :param invalidates_cache: when set to `True`, a process exiting
            with this exit code will not be considered for caching
        """
        if invalidates_cache is None:
            invalidates_cache = (isinstance(status, int) and status < 400)

        if 'invalidates_cache' in ExitCode._fields:
            super(CalcJobProcessSpec,
                  self).exit_code(status=status, label=label, message=message, invalidates_cache=invalidates_cache)
        else:
            super(CalcJobProcessSpec, self).exit_code(status=status, message=message, label=label)


class CalcJob(_BaseCalcJob):  # pylint: disable=abstract-method
    """Custom CalcJob class for aiida-quantumespresso calculations."""

    _spec_class = CalcJobProcessSpec
