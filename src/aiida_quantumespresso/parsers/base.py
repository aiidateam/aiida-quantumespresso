# -*- coding: utf-8 -*-
"""Defines a `Parser` base class for `aiida-quantumespresso`.

All `Parser` implementations in `aiida-quantumespresso` must use this base class, not `aiida.parsers.Parser`.
"""
from aiida.parsers import Parser as _BaseParser

__all__ = ('Parser',)


class Parser(_BaseParser):  # pylint: disable=abstract-method
    """Custom `Parser` class for `aiida-quantumespresso` parser implementations."""

    def emit_logs(self, logging_dictionaries, ignore=None):
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

    def exit(self, exit_code):
        """Log the exit message of the give exit code with level `ERROR` and return the exit code.

        This is a utility function if one wants to return from the parse method and automically add the exit message
        associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

        :param exit_code: an `ExitCode`
        :return: the exit code
        """
        self.logger.error(exit_code.message)
        return exit_code
