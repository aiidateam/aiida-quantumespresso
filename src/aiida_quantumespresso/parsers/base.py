# -*- coding: utf-8 -*-
"""Defines a `Parser` base class for `aiida-quantumespresso`.

All `Parser` implementations in `aiida-quantumespresso` must use this base class, not `aiida.parsers.Parser`.
"""
import abc
import re
import typing

from aiida.common import AttributeDict
from aiida.engine import ExitCode
from aiida.parsers import Parser

from aiida_quantumespresso.parsers.parse_raw.base import convert_qe_time_to_sec
from aiida_quantumespresso.utils.mapping import get_logging_container

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

    def _parse_stdout_from_retrieved(self, **kwargs) -> typing.Tuple[str, dict, AttributeDict]:
        """Retrieve and parse the ``stdout`` content of a Quantum ESPRESSO calculation.

        :returns: size 3 tuple with the stdout content, parsed data and log messages
        """
        logs = get_logging_container()

        filename_stdout = self.node.get_option('output_filename')

        if filename_stdout not in self.retrieved.base.repository.list_object_names():
            logs.error.append('ERROR_OUTPUT_STDOUT_MISSING')
            return {}, logs

        try:
            with self.retrieved.open(filename_stdout, 'r') as handle:
                stdout = handle.read()
        except OSError:
            logs.error.append('ERROR_OUTPUT_STDOUT_READ')
            return {}, logs

        try:
            parsed_data, stdout_logs = self.parse_stdout(stdout, **kwargs)
        except Exception as exception:
            logs.error.append('ERROR_OUTPUT_STDOUT_PARSE')
            logs.error.append(exception)
            return {}, logs

        for log_level, log_items in stdout_logs.items():
            logs[log_level].extend(log_items)

        return parsed_data, logs

    @classmethod
    def parse_stdout(cls, stdout: str) -> typing.Tuple[dict, AttributeDict]:
        """Parse the ``stdout`` content of a Quantum ESPRESSO calculation.

        This function only checks for basic content like JOB DONE, errors with %%%%% etc.

        :param stdout: the stdout content as a string.
        :returns: tuple of two dictionaries, with the parsed data and log messages, respectively.
        """
        logs = get_logging_container()
        parsed_data = {}

        if not re.search(r'JOB DONE', stdout):
            logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

        code_match = re.search(r'Program\s(?P<code_name>[A-Z|\_|\d]+)\s(?P<code_version>v\.[\d\.|a-z|A-Z]+)\s', stdout)

        if code_match:

            code_name = code_match.groupdict()['code_name']
            parsed_data['code_version'] = code_match.groupdict()['code_version']

            wall_match = re.search(fr'{code_name}\s+:[\s\S]+\s+(?P<wall_time>[.\d|s|m|d|h]+)\sWALL', stdout)

            if wall_match:
                parsed_data['wall_time'] = wall_match.groupdict()['wall_time']

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

    def _emit_logs(self, logging_dictionaries: AttributeDict, ignore: list = None) -> None:
        """Emit the messages in one or multiple "log dictionaries" through the logger of the parser.

        A log dictionary is expected to have the following structure: each key must correspond to a log level of the
        python logging module, e.g. `error` or `warning` and its values must be a list of string messages. The method
        will loop over all log dictionaries and emit the messages it contains with the log level indicated by the key.

        Example log dictionary structure::

            logs = {
                'warning': ['Could not parse the `etot_threshold` variable from the stdout.'],
                'error': ['Self-consistency was not achieved']
            }

        :param logging_dictionaries: log dictionaries
        :param ignore: list of log messages to ignore
        """
        ignore = ignore or []

        if not isinstance(logging_dictionaries, (list, tuple)):
            logging_dictionaries = [logging_dictionaries]

        for logs in logging_dictionaries:
            for level, messages in logs.items():
                for message in messages:

                    if message is None:
                        continue

                    stripped = message.strip()

                    if not stripped or stripped in ignore:
                        continue

                    try:
                        getattr(self.logger, level)(stripped)
                    except AttributeError:
                        pass

    def _exit(self, exit_code: ExitCode) -> ExitCode:
        """Log the exit message of the give exit code with level `ERROR` and return the exit code.

        This is a utility function if one wants to return from the parse method and automically add the exit message
        associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

        :param exit_code: an `ExitCode`
        :return: the exit code
        """
        self.logger.error(exit_code.message)
        return exit_code
