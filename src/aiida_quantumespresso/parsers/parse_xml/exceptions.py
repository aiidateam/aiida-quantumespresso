# -*- coding: utf-8 -*-
"""Exceptions for the XML parsing module of Quantum ESPRESSO."""


class XMLParseError(Exception):
    """Raised when the XML output could not be parsed."""


class XMLUnsupportedFormatError(Exception):
    """Raised when the XML output has an unsupported format."""
