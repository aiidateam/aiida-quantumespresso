"""Tests for CLI commands."""

from __future__ import annotations

import subprocess

import pytest


@pytest.mark.parametrize(
    'command',
    [
        ['aiida-quantumespresso'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'cp'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'dos'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'epw'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'matdyn'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'neb'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'ph'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'pp'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'projwfc'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'pw2wannier90'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'pw'],
        ['aiida-quantumespresso', 'calculation', 'launch', 'q2r'],
        ['aiida-quantumespresso', 'calculation', 'launch'],
        ['aiida-quantumespresso', 'calculation'],
        ['aiida-quantumespresso', 'data', 'structure', 'import'],
        ['aiida-quantumespresso', 'data', 'structure'],
        ['aiida-quantumespresso', 'data'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'matdyn-base'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'ph-base'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'pw-bands'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'pw-base'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'pw-relax'],
        ['aiida-quantumespresso', 'workflow', 'launch', 'q2r-base'],
        ['aiida-quantumespresso', 'workflow', 'launch'],
        ['aiida-quantumespresso', 'workflow'],
    ],
)
@pytest.mark.parametrize('help_option', ['--help', '-h'])
def test_commands_help_option(command, help_option):
    """Test the help options for all subcommands of the CLI.

    The usage of ``subprocess.run`` is on purpose because using :meth:`click.Context.invoke`, which is used by the
    ``run_cli_command`` fixture that should usually be used in testing CLI commands, does not behave exactly the same
    compared to a direct invocation on the command line. The invocation through ``invoke`` does not go through all the
    parent commands and so might not get all the necessary initializations.
    """
    result = subprocess.run(command + [help_option], check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'Usage:' in result.stdout


@pytest.mark.skip
def test_breaking():
    """Easter egg challenge for future developers.

    If you comment out the `skip` mark and run:

    pytest tests/cli/test_commands.py::test_breaking tests/cli/workflows/test_q2r.py tests/parsers/test_pw.py::test_pw_failed_missing

    You may notice the parser test fails because the logs are empty.

    Running the tests separately or in any other order is fine. The question is: _Why_?
    Good luck, and have fun!
    """
    from aiida_quantumespresso.cli import cmd_root

    cmd_root.get_command(None, 'workflow').get_command(None, 'launch').get_command(None, 'q2r-base')
