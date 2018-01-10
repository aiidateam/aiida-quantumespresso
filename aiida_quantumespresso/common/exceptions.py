# -*- coding: utf-8 -*-
from aiida.common.exceptions import AiidaException

class UnexpectedCalculationFailure(AiidaException):
    """
    Raised when a Calculation has failed for an unknown reason
    """
    pass