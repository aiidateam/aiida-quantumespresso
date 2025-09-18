# -*- coding: utf-8 -*-
"""Command line interface commands for setting up `code`s for Quantum ESPRESSO executables."""

from aiida import orm
from aiida.cmdline.params import arguments
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent
import click

from aiida_quantumespresso.tools.code_setup import get_code_label, get_executable_paths


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
@click.option(
    '--prepend-text',
    default='',
    required=False,
    help=(
        'Bash commands to prepend to the executable call in all submit scripts '
        'for the codes. Pass a string directly, or use the `--interactive` '
        'option to open an editor.'
    ),
)
@click.option(
    '--append-text',
    default='',
    required=False,
    help=(
        'Bash commands to append to the executable call in all submit scripts '
        'for the codes. Pass a string directly, or use the `--interactive` '
        'option to open an editor.'
    ),
)
@click.option(
    '--interactive',
    '-I',
    is_flag=True,
    help='Open an editor to edit the prepend and append text.',
)
def create_codes_cmd(computer, executables, directory, label_template, prepend_text, append_text, interactive):
    """Automatically create an `orm.Code` instance for Quantum ESPRESSO executables.

    Specify the target `orm.Computer` and a single executable or a list of Quantum ESPRESSO executables to
    create codes for. You can provide multiple executables separated by a
    space, e.g. `pw.x dos.x`.

    Use `--directory`/`-d` to point to the directory containing the executables.
    If not provided, the command will try to find the executables in
    the `PATH` on the (remote) computer. Consider adding a `--prepend-text` if modules
    need to be loaded to find the executables. This can also be done with an editor via
    the `--interactive` option.
    """
    if interactive:
        prepend_text = click.edit(prepend_text or '# Enter PREPEND text here...') or prepend_text
        append_text = click.edit(append_text or '# Enter APPEND text here...') or append_text

    user = orm.User.collection.get_default()

    try:
        computer.get_authinfo(user)
    except NotExistent:
        echo.echo_critical(f'Computer<{computer.label}> is not yet configured for user<{user.email}>')

    try:
        executable_path_mapping = get_executable_paths(executables, computer, prepend_text, directory)
    except (FileNotFoundError, ValueError) as exc:
        echo.echo_critical(exc)

    for executable, exec_path in executable_path_mapping.items():
        existing_label, label = get_code_label(label_template=label_template, executable=executable, computer=computer)
        if existing_label:
            echo.echo_warning(f'Code with label<{existing_label}> already exists on Computer<{computer.label}>.')
            if not click.confirm(f'Do you want to add another instance with label {label}?'):
                continue

        code = orm.InstalledCode(
            label=label,
            computer=computer,
            filepath_executable=exec_path,
            default_calc_job_plugin=f'quantumespresso.{executable.split(".")[0]}',
            prepend_text=prepend_text,
        )
        code.store()
        echo.echo_success(
            f'Code<{code.label}> for {executable} created with pk<{code.pk}> on Computer<{computer.label}>'
        )
