"""Command for setting up codes for Quantum ESPRESSO executables."""

import click
from aiida import orm
from aiida.cmdline.params import arguments
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent

from aiida_quantumespresso.tools.setup import get_code_label, get_executable_paths, VALID_QE_EXECUTABLES

PREPEND_APPEND_TEMPLATE = (
    '#==========================================================================\n'
    '#= Enter your {} text below. Lines starting with "#=" will be ignored.\n'
    '#=========================================================================='
)


@arguments.COMPUTER()
@click.argument(
    'executables',
    nargs=-1,
    required=False,
)
@click.option(
    '--directory',
    '-d',
    help='Absolute path to directory where the executable(s) is (are) located.',
)
@click.option(
    '--label-template',
    '-l',
    help=(
        'Label template for the code instance. Use curly brackets to reference '
        'the executable label, e.g. `qe-{}` will create a `Code` with label '
        '`qe-pw` for `pw.x`. Defaults to the executable name without `.x` suffix., e.g. `pw` for `pw.x`.'
    ),
    default='',
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
    '-i',
    is_flag=True,
    help='Open an editor to edit the prepend and append text.',
)
@click.option(
    '--all',
    '-a',
    'setup_all',
    is_flag=True,
    help='Set up codes for all supported Quantum ESPRESSO executables.',
)
@click.option(
    '--ignore-missing',
    is_flag=True,
    help='Ignore missing executables and do not raise an error.',
)
@click.option(
    '--skip-existing',
    '-s',
    is_flag=True,
    help='Skip creating codes if they already exist on the computer.',
)
def setup_codes_cmd(
    computer,
    executables,
    directory,
    label_template,
    prepend_text,
    append_text,
    interactive,
    setup_all,
    ignore_missing,
    skip_existing,
):
    """Set up codes for Quantum ESPRESSO executables.

    Specify the target `orm.Computer` and a single executable or a list of Quantum ESPRESSO executables to
    create codes for. You can provide multiple executables separated by a
    space, e.g. `pw.x dos.x`. To set up codes for all supported Quantum ESPRESSO executables, use the `--all`/`-a` option.

    Use `--directory`/`-d` to point to the directory containing the executables.
    If not provided, the command will try to find the executables in
    the `PATH` on the (remote) computer. Consider adding a `--prepend-text` if modules
    need to be loaded to find the executables. This can also be done with an editor via
    the `--interactive` option.
    """
    if bool(executables) == setup_all:
        if setup_all:
            echo.echo_critical('Please provide either executables or use the `--all/-a` option, not both.')
        echo.echo_critical('Please provide at least one executable or use the `--all/-a` option to set up all codes.')

    if interactive:
        prepend_text = click.edit(prepend_text or PREPEND_APPEND_TEMPLATE.format('PREPEND')) or prepend_text
        prepend_text = '\n'.join([line for line in prepend_text.splitlines() if not line.strip().startswith('#=')])
        append_text = click.edit(append_text or PREPEND_APPEND_TEMPLATE.format('APPEND')) or append_text
        append_text = '\n'.join([line for line in append_text.splitlines() if not line.strip().startswith('#=')])

    user = orm.User.collection.get_default()

    try:
        computer.get_authinfo(user)
    except NotExistent:
        echo.echo_critical(f'Computer<{computer.label}> is not yet configured for user<{user.email}>')

    if setup_all:
        ignore_missing = True
        executables = VALID_QE_EXECUTABLES

    try:
        executable_path_mapping = get_executable_paths(executables, computer, prepend_text, directory, ignore_missing)
    except (FileNotFoundError, ValueError) as exc:
        echo.echo_critical(exc)

    if ignore_missing:
        missing = False
        for exe in executables:
            if exe not in executable_path_mapping:
                echo.echo_warning(
                    f'Executable<{exe}> not found on Computer<{computer.label}>. Skipping setup for this executable.'
                )
                missing = True
        if missing:
            echo.echo('')

    for executable, exec_path in executable_path_mapping.items():
        existing_label, label = get_code_label(label_template=label_template, executable=executable, computer=computer)
        if existing_label:
            if skip_existing:
                echo.echo_warning(
                    f'Code<{existing_label}> for {executable} already exists on Computer<{computer.label}>. '
                    'Skipping setup for this executable.'
                )
                continue
            echo.echo_warning(f'Code with label<{existing_label}> already exists on Computer<{computer.label}>.')
            if not click.confirm(f'Do you want to add another instance with label {label}?'):
                continue

        code = orm.InstalledCode(
            label=label,
            computer=computer,
            filepath_executable=exec_path,
            default_calc_job_plugin=f'quantumespresso.{executable.split(".")[0]}',
            prepend_text=prepend_text,
            append_text=append_text,
        )
        code.store()
        echo.echo_success(
            f'Code<{code.label}> for {executable} created with pk<{code.pk}> on Computer<{computer.label}>'
        )
