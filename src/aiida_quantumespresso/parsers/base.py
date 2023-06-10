# -*- coding: utf-8 -*-
"""Defines a `Parser` base class for `aiida-quantumespresso`.

All `Parser` implementations in `aiida-quantumespresso` must use this base class, not `aiida.parsers.Parser`.
"""
from __future__ import annotations

import abc
import re
from typing import Optional, Tuple

from aiida.common import AttributeDict
from aiida.engine import ExitCode
from aiida.parsers import Parser

from aiida_quantumespresso.parsers.parse_raw.base import convert_qe_time_to_sec

__all__ = ('BaseParser',)


class BaseParser(Parser, metaclass=abc.ABCMeta):
    """Custom ``Parser`` class for ``aiida-quantumespresso`` parser implementations."""

    class_error_map = {}
    class_warning_map = {}

    base_error_map = {
        'Maximum CPU time exceeded': 'ERROR_OUT_OF_WALLTIME',
    }
    base_warning_map = {
        'Warning:': None,
        'DEPRECATED:': None,
    }
    success_string = 'JOB DONE'

    @classmethod
    def get_error_map(cls):
        """The full error map of the parser class."""
        error_map = cls.base_error_map.copy()
        error_map.update(cls.class_error_map)
        return error_map

    @classmethod
    def get_warning_map(cls):
        """The full warning map of the parser class."""
        warning_map = cls.base_warning_map.copy()
        warning_map.update(cls.class_warning_map)
        return warning_map

    def parse_stdout_from_retrieved(self, logs: AttributeDict) -> Tuple[str, dict, AttributeDict]:
        """Read and parse the ``stdout`` content of a Quantum ESPRESSO calculation.

        :param logs: Logging container that will be updated during parsing.
        :returns: size 3 tuple: (``stdout`` content, parsed data, updated logs).
        """
        filename_stdout = self.node.get_option('output_filename')

        if filename_stdout not in self.retrieved.base.repository.list_object_names():
            logs.error.append('ERROR_OUTPUT_STDOUT_MISSING')
            return '', {}, logs

        try:
            with self.retrieved.open(filename_stdout, 'r') as handle:
                stdout = handle.read()
        except OSError as exception:
            logs.error.append('ERROR_OUTPUT_STDOUT_READ')
            logs.error.append(exception)
            return '', {}, logs

        try:
            parsed_data, logs = self._parse_stdout_base(stdout, logs)
        except Exception as exception:
            logs.error.append('ERROR_OUTPUT_STDOUT_PARSE')
            logs.error.append(exception)
            return stdout, {}, logs

        return stdout, parsed_data, logs

    def emit_logs(self, logs: list[AttributeDict] | tuple[AttributeDict] | AttributeDict, ignore: list = None) -> None:
        """Emit the messages in one or multiple "log dictionaries" through the logger of the parser.

        A log dictionary is expected to have the following structure: each key must correspond to a log level of the
        python logging module, e.g. `error` or `warning` and its values must be a list of string messages. The method
        will loop over all log dictionaries and emit the messages it contains with the log level indicated by the key.

        Example log dictionary structure::

            logs = {
                'warning': ['Could not parse the `etot_threshold` variable from the stdout.'],
                'error': ['Self-consistency was not achieved']
            }

        :param logs: log dictionaries
        :param ignore: list of log messages to ignore
        """
        ignore = ignore or []

        if not isinstance(logs, (list, tuple)):
            logs = [logs]

        for logs in logs:
            for level, messages in logs.items():
                for message in messages:

                    stripped = message.strip()

                    if stripped in ignore:
                        continue

                    getattr(self.logger, level)(stripped)

    def check_base_errors(self, logs: AttributeDict) -> Optional[ExitCode]:
        """Check the ``logs`` for the following "basic" parsing error and return a (formatted) version:

        * ``ERROR_OUTPUT_STDOUT_MISSING``
        * ``ERROR_OUTPUT_STDOUT_READ``
        * ``ERROR_OUTPUT_STDOUT_PARSE``

        These errors mean that there is no ``stdout`` to parse.

        The ``ERROR_OUTPUT_STDOUT_INCOMPLETE`` error is not checked here because in this case there might still be
        useful information in the ``stdout``.
        """

        for exit_code in [
            'ERROR_OUTPUT_STDOUT_MISSING',
        ]:
            if exit_code in logs.error:
                return self.exit_codes.get(exit_code)

        # These exit codes have additional information that needs to be formatted in the message.
        for exit_code in [
            'ERROR_OUTPUT_STDOUT_READ',
            'ERROR_OUTPUT_STDOUT_PARSE'
        ]:
            if exit_code in logs.error:
                exception = logs.error[logs.index(exit_code) + 1]
                return self.exit_codes.get(exit_code).format(exception=exception)

    def exit(self, exit_code: ExitCode | None = None, logs: AttributeDict | None = None) -> ExitCode:
        """Log all messages in the ``logs`` as well as the ``exit_code`` message and return the correct exit code.

        This is a utility function if one wants to return from the parse method and automically add the ``logs`` and
        exit message associated to and exit code as a log message to the node: e.g.
        ``return self._exit(self.exit_codes.LABEL))``

        If no ``exit_code`` is provided, the method will check if an ``exit_status`` has already been set on the node
        and return the corresponding ``ExitCode`` in this case. If not, ``ExitCode(0)`` is returned.

        :param logs: log dictionaries
        :param exit_code: an ``ExitCode``
        :return: The correct exit code
        """
        if logs:
            self.emit_logs(logs)

        if exit_code is not None:
            self.logger.error(exit_code.message)
        elif self.node.exit_status is not None:
            exit_code = ExitCode(self.node.exit_status, self.node.exit_message)
        else:
            exit_code = ExitCode(0)

        return exit_code

    @classmethod
    def _parse_stdout_base(cls, stdout: str, logs: AttributeDict) -> Tuple[dict, AttributeDict]:
        """Parse the ``stdout`` content of a Quantum ESPRESSO calculation.

        This function only checks for basic content like JOB DONE, errors with %%%%% etc, but can be overridden to
        parse more data from the ``stdout``.

        :param stdout: the stdout content as a string.
        :returns: tuple of two dictionaries, with the parsed data and log messages, respectively.
        """
        parsed_data = {}

        if not re.search(cls.success_string, stdout):
            logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

        code_match = re.search(
            r'Program\s(?P<code_name>[A-Z|a-z|\_|\d]+)\sv\.(?P<code_version>[\d\.|a-z|A-Z]+)\s', stdout
        )
        if code_match:
            code_name = code_match.groupdict()['code_name']
            parsed_data['code_version'] = code_match.groupdict()['code_version']

            wall_match = re.search(fr'{code_name}\s+:[\s\S]+CPU\s+(?P<wall_time>[\s.\d|s|m|d|h]+)\sWALL', stdout)

            if wall_match:
                try:
                    parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(wall_match.groupdict()['wall_time'])
                except ValueError:
                    logs.warnings.append('Unable to convert wall time from `stdout` to seconds.')

        # Look for typical Quantum ESPRESSO error messages between %%%%%-lines that are not in our error map
        if re.search(r'\%\%\%\%\%', stdout):  # Note: using e.g. `\%{5}` is significantly slower
            for error_message in set(re.split(r'\%\%\%\%\%\n', stdout)[1::2]):

                if not any(error_marker in error_message for error_marker in cls.get_error_map().keys()):
                    logs.error.append(error_message.rstrip('\n%'))

        # Look for error messages in general
        for error_marker, error, in cls.get_error_map().items():
            if re.search(fr'{error_marker}', stdout):
                logs.error.append(error)

        # Look for lines with warnings from the `warning_map`
        for warning_marker, warning in cls.get_warning_map().items():
            for warning_message in set(re.findall(fr'({warning_marker}.+)\n', stdout)):
                if warning is not None:
                    logs.warning.append(warning)
                else:
                    logs.warning.append(warning_message)

        return parsed_data, logs
