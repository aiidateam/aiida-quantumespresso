# -*- coding: utf-8 -*-
"""Tests for CLI commands."""
from click import Context, Group
from aiida_quantumespresso.cli import cmd_root


def test_commands():
    """Test that all commands in ``cmd_root`` are reachable and can print the help message.

    This doesn't guarantee that the command works but at least that it can be successfully called and there are no
    import errors or other basic problems.
    """

    def recursively_print_help(ctx):
        assert isinstance(ctx.get_help(), str)

        if isinstance(ctx.command, Group):
            for subcommand in ctx.command.commands.values():
                ctx.command = subcommand
                recursively_print_help(ctx)

    ctx = Context(cmd_root)
    recursively_print_help(ctx)
