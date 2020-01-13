# -*- coding: utf-8 -*-

"""
Defines custom base classes for CalcJob and WorkChain which automatically
set the 'invalidates_cache' attribute of exit codes based on the status
integer.
"""

from aiida.engine import ExitCode
from aiida.engine import CalcJob as BaseCalcJob
from aiida.engine import WorkChain as BaseWorkChain
from aiida.engine.processes.process_spec import ProcessSpec as BaseProcessSpec


__all__ = ('ProcessSpec', 'CalcJob', 'WorkChain')

class ProcessSpec(BaseProcessSpec):

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
            super(ProcessSpec,
                  self).exit_code(status=status, label=label, message=message, invalidates_cache=invalidates_cache)
        else:
            super(ProcessSpec, self).exit_code(status=status, message=message, label=label)


class CalcJob(BaseCalcJob):
    _spec_class = ProcessSpec


class WorkChain(BaseWorkChain):
    _spec_class = ProcessSpec
