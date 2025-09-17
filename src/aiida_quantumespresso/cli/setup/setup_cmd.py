# -*- coding: utf-8 -*-
"""Command line interface commands for setting up `code`s for Quantum ESPRESSO executables."""

from aiida.cmdline.params import arguments
from aiida.cmdline.params.options.commands import code as options_code
import click

from aiida_quantumespresso.tools.code_setup import _create_codes


@arguments.COMPUTER()
@click.argument(
    'executables',
    nargs=-1,
    required=True,
)
@click.option('--directory', '-d', help='Absolute path to directory where the executable(s) is (are) located.')
@click.option(
    '--label-template',
    '-l',
    help=(
        'Label template for the code instance. Use curly brackets to reference '
        'the executable label, e.g. `qe-{}` will create a `Code` with label '
        '`qe-pw` for `pw.x`. Defaults to the executable name without `.x` suffix., e.g. `pw` for `pw.x`.'
    ),
    default=''
)
@options_code.PREPEND_TEXT()
@options_code.APPEND_TEXT()
def create_codes_cmd(computer, executables, directory, label_template, **kwargs):
    """Automatically create an `orm.Code` instance for Quantum ESPRESSO executables.

    Specify the target `orm.Computer` and a single executable or a list of Quantum ESPRESSO executables to
    create codes for. You can provide multiple executables separated by a
    space, e.g. `pw.x dos.x`.

    Specify a directory where the executables are located with `--directory` or
    `-d`. If not provided, the command will try to find the executables in
    the `PATH` on the (remote) computer. Consider adding a prepend text if modules
    need to be loaded to find the executables.
    """
    prepend_text = kwargs.pop('prepend_text', None)
    append_text = kwargs.pop('append_text', None)

    _create_codes(
        computer=computer,
        executables=executables,
        directory=directory,
        label_template=label_template,
        prepend_text=prepend_text,
        append_text=append_text,
        on_conflict='prompt',
    )
