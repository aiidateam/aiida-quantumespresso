# -*- coding: utf-8 -*-
"""Utilities for `WorkChain` implementations."""
from __future__ import absolute_import

from collections import namedtuple
from functools import wraps

from aiida.engine import ExitCode

ErrorHandler = namedtuple('ErrorHandler', 'priority method')
"""
A namedtuple to define an error handler for a :class:`~aiida.engine.processes.workchains.workchain.WorkChain`.

The priority determines in which order the error handling methods are executed, with
the higher priority being executed first. The method defines an unbound WorkChain method
that takes an instance of a :class:`~aiida.orm.nodes.process.calculations.calcjob.CalcJobNode`
as its sole argument. If the condition of the error handler is met, it should return an :class:`.ErrorHandlerReport`.

:param priority: integer denoting the error handlers priority
:param method: the workchain class method
"""

ErrorHandlerReport = namedtuple('ErrorHandlerReport', 'is_handled do_break exit_code')
ErrorHandlerReport.__new__.__defaults__ = (False, False, ExitCode())
"""
A namedtuple to define an error handler report for a :class:`~aiida.engine.processes.workchains.workchain.WorkChain`.

This namedtuple should be returned by an error handling method of a workchain instance if
the condition of the error handling was met by the failure mode of the calculation.
If the error was appriopriately handled, the 'is_handled' field should be set to `True`,
and `False` otherwise. If no further error handling should be performed after this method
the 'do_break' field should be set to `True`

:param is_handled: boolean, set to `True` when an error was handled, default is `False`
:param do_break: boolean, set to `True` if no further error handling should be performed, default is `False`
:param exit_code: an instance of the :class:`~aiida.engine.processes.exit_code.ExitCode` tuple
"""


def register_error_handler(cls, priority):
    """
    Decorator that will turn any function in an error handler for workchain that inherits from
    the :class:`.BaseRestartWorkChain`. The function expects two arguments, a workchain class and a priortity.
    The decorator will add the function as a class method to the workchain class and add an :class:`.ErrorHandler`
    tuple to the :attr:`.BaseRestartWorkChain._error_handlers` attribute of the workchain. During failed calculation
    handling the :meth:`.inspect_calculation` outline method will call the `_handle_calculation_failure` which will loop
    over all error handler in the :attr:`.BaseRestartWorkChain._error_handlers`, sorted with respect to the priority in
    reverse. If the workchain class defines a :attr:`.BaseRestartWorkChain._verbose` attribute and is set to `True`, a
    report message will be fired when the error handler is executed.

    Requirements on the function signature of error handling functions. The function to which the
    decorator is applied needs to take two arguments:

        * `self`: This is the instance of the workchain itself
        * `calculation`: This is the calculation that failed and needs to be investigated

    The function body should usually consist of a single conditional that checks the calculation if
    the error that it is designed to handle is applicable. Although not required, it is advised that
    the function return an :class:`.ErrorHandlerReport` tuple when its conditional was met. If an error was handled
    it should set `is_handled` to `True`. If no other error handlers should be considered set `do_break` to `True`.

    :param cls: the workchain class to register the error handler with
    :param priority: an integer that defines the order in which registered handlers will be called
        during the handling of a failed calculation. Higher priorities will be handled first
    """

    def error_handler_decorator(handler):
        """Decorator to dynamically register an error handler to a `WorkChain` class."""

        @wraps(handler)
        def error_handler(self, calculation):
            """Wrapped error handler to add a log to the report if the handler is called and verbosity is turned on."""
            if hasattr(cls, '_verbose') and cls._verbose:  # pylint: disable=protected-access
                self.report('({}){}'.format(priority, handler.__name__))
            return handler(self, calculation)

        setattr(cls, handler.__name__, error_handler)

        if not hasattr(cls, '_error_handlers'):
            cls._error_handlers = []  # pylint: disable=protected-access
        cls._error_handlers.append(ErrorHandler(priority, error_handler))  # pylint: disable=protected-access

        return error_handler

    return error_handler_decorator
