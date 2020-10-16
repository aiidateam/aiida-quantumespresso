# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Fixtures for the command line interface."""
import pytest


def mock_launch_process(*_, **__):
    """Mock the :meth:`~aiida_quantumespresso.cli.utilslaunch.launch_process` to be a no-op."""
    return


@pytest.fixture
def run_cli_command():
    """Run a `click` command with the given options.

    The call will raise if the command triggered an exception or the exit code returned is non-zero.
    """

    def _run_cli_command(command, options=None, raises=None):
        """Run the command and check the result.

        :param command: the command to invoke
        :param options: the list of command line options to pass to the command invocation
        :param raises: optionally an exception class that is expected to be raised
        """
        import traceback
        from click.testing import CliRunner

        runner = CliRunner()
        result = runner.invoke(command, options or [])

        if raises is not None:
            assert result.exception is not None, result.output
            assert result.exit_code != 0
        else:
            assert result.exception is None, ''.join(traceback.format_exception(*result.exc_info))
            assert result.exit_code == 0, result.output

        result.output_lines = [line.strip() for line in result.output.split('\n') if line.strip()]

        return result

    return _run_cli_command


@pytest.fixture
def run_cli_process_launch_command(run_cli_command, monkeypatch):
    """Run a process launch command with the given options.

    The call will raise if the command triggered an exception or the exit code returned is non-zero.

    :param command: the command to invoke
    :param options: the list of command line options to pass to the command invocation
    :param raises: optionally an exception class that is expected to be raised
    """

    def _inner(command, options=None, raises=None):
        """Run the command and check the result."""
        from aiida_quantumespresso.cli.utils import launch
        monkeypatch.setattr(launch, 'launch_process', mock_launch_process)
        return run_cli_command(command, options, raises)

    return _inner
