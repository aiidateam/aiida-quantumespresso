# -*- coding: utf-8 -*-
# pylint: disable=wrong-import-position
"""Module for the command line interface."""
from aiida.cmdline.params import options, types
import click


class VerbosityGroup(click.Group):
    """Custom command group that automatically adds the ``VERBOSITY`` option to all subcommands."""

    @staticmethod
    def add_verbosity_option(cmd):
        """Apply the ``verbosity`` option to the command, which is common to all ``verdi`` commands."""
        if 'verbosity' not in [param.name for param in cmd.params]:
            cmd = options.VERBOSITY()(cmd)

        return cmd

    def group(self, *args, **kwargs):
        """Ensure that sub command groups use the same class but do not override an explicitly set value."""
        kwargs.setdefault('cls', self.__class__)
        return super().group(*args, **kwargs)

    def get_command(self, ctx, cmd_name):
        """Return the command that corresponds to the requested ``cmd_name``.

        This method is overridden from the base class in order to automatically add the verbosity option.

        Note that if the command is not found and ``resilient_parsing`` is set to True on the context, then the latter
        feature is disabled because most likely we are operating in tab-completion mode.
        """
        cmd = super().get_command(ctx, cmd_name)

        if cmd is not None:
            return self.add_verbosity_option(cmd)

        if ctx.resilient_parsing:
            return None

        return ctx.fail(f'`{cmd_name}` is not a {self.name} command.')


@click.group('aiida-quantumespresso', context_settings={'help_option_names': ['-h', '--help']})
@options.PROFILE(type=types.ProfileParamType(load_profile=True), expose_value=False)
@options.VERBOSITY()
def cmd_root():
    """CLI for the `aiida-quantumespresso` plugin."""


from .calculations import cmd_calculation
from .data import cmd_structure
from .workflows import cmd_workflow
