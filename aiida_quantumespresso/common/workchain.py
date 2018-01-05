# -*- coding: utf-8 -*-
from collections import namedtuple

"""
A namedtuple to define an error handler for a WorkChain. It defines two fields:

    * priority: an integer
    * method: a WorkChain method

The priority determines in which order the error handling methods are executed, with
the higher priority being executed first. The method defines an unbound WorkChain method
that takes an instance of a Calculation as its sole argument. If the condition of the
error handler is met, it should return an ErrorHandlerReport
"""
ErrorHandler = namedtuple('ErrorHandler', 'priority method')


"""
A namedtuple to define an error handler report for a WorkChain. It defines two fields:

    * is_handled: a boolean
    * do_break: a boolean

This namedtuple should be returned by an error handling method of WorkChain class if
the condition of the error handling was met by the failure mode of the Calculation.
If the error was appriopriately handled, the 'is_handled' field should be set to True,
and False otherwise. If no further error handling should be performed after this method
the 'do_break' field should be set to True
"""
ErrorHandlerReport = namedtuple('ErrorHandlerReport', 'is_handled do_break')