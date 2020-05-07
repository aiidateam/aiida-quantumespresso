# -*- coding: utf-8 -*-
"""Exceptions specific to `aiida-quantumespresso`."""
from aiida.common import AiidaException


class UnexpectedCalculationFailure(AiidaException):
    """Raised when a calculation job has failed for an unexpected or unrecognized reason."""
