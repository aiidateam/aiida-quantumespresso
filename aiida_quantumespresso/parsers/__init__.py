# -*- coding: utf-8 -*-
from aiida.common import OutputParsingError


class QEOutputParsingError(OutputParsingError):
    """Exception raised when there is a parsing error in the QE parser."""
    pass
