# -*- coding: utf-8 -*-
"""Defines CalcJob and WorkChain base classes for aiida-quantumespresso.

The custom CalcJob and WorkChain base classes automatically set the
'invalidates_cache' attribute of exit codes based on the status integer.
"""

from __future__ import absolute_import

from aiida.engine import ExitCode
from aiida.engine import CalcJob
from aiida.engine import WorkChain
from aiida.engine.processes.process_spec import CalcJobProcessSpec
from aiida.engine.processes.workchains.workchain import WorkChainSpec

__all__ = ('QuantumEspressoCalcJob', 'QuantumEspressoWorkChain')


class InvalidatesCacheProcessSpecMixin(object):  # pylint: disable=too-few-public-methods
    """Implements automatically setting the `invalidates_cache` flag.

    Mix-in class which adds the feature of automatically setting the
    `invalidates_cache` flag to `True` for an exit status smaller than
    400, unless explicitly overriden.
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
            super(InvalidatesCacheProcessSpecMixin,
                  self).exit_code(status=status, label=label, message=message, invalidates_cache=invalidates_cache)
        else:
            super(InvalidatesCacheProcessSpecMixin, self).exit_code(status=status, message=message, label=label)


class QuantumEspressoCalcJobProcessSpec(InvalidatesCacheProcessSpecMixin, CalcJobProcessSpec):  # pylint: disable=missing-docstring
    pass


class QuantumEspressoWorkChainSpec(InvalidatesCacheProcessSpecMixin, WorkChainSpec):  # pylint: disable=missing-docstring
    pass


class QuantumEspressoCalcJob(CalcJob):  # pylint: disable=missing-docstring,abstract-method
    _spec_class = QuantumEspressoCalcJobProcessSpec


class QuantumEspressoWorkChain(WorkChain):  # pylint: disable=missing-docstring
    _spec_class = QuantumEspressoWorkChainSpec
