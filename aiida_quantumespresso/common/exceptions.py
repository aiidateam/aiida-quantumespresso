# -*- coding: utf-8 -*-
"""Exceptions specific to `aiida-quantumespresso`."""
from __future__ import absolute_import

from aiida.common import AiidaException


class UnexpectedCalculationFailure(AiidaException):
    """Raised when a calculation job has failed for an unexpected or unrecognized reason."""
