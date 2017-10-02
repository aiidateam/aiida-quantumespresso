# -*- coding: utf-8 -*-
from aiida.common.exceptions import AiidaException

class UnexpectedFailure(AiidaException):
    """
    Raised when a PwCalculation has failed for an unknown reason
    """
    pass