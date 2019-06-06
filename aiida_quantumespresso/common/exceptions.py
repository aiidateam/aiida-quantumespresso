# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.common import AiidaException


class UnexpectedCalculationFailure(AiidaException):
    """Raised when a calculation has failed for an unknown reason."""
    pass
