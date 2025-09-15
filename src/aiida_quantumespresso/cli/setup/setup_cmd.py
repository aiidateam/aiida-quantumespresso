# -*- coding: utf-8 -*-
"""Command line interface commands for setting up `code`s for Quantum ESPRESSO executables."""
import pathlib
import re
from typing import List, Tuple, Union

from aiida import orm
from aiida.cmdline.params import arguments
from aiida.cmdline.params.options.commands import code as options_code
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent
import click


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
def create_code(computer, executables, directory, label_template, **kwargs):
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

    executable_paths = _get_executable_paths(
        prepend_text=prepend_text, executables=executables, computer=computer, directory=directory, **kwargs
    )

    for executable, exec_path in executable_paths:
        existing_label, label = _get_code_label(label_template=label_template, executable=executable, computer=computer)

        if existing_label:
            echo.echo_warning(f'Code with label<{existing_label}> already '
                              f'exists on Computer<{computer.label}>.')
            if not click.confirm(f'Do you want to add another instance with label {label}?'):
                continue

        code = orm.InstalledCode(
            label=label,
            computer=computer,
            filepath_executable=exec_path,
            default_calc_job_plugin=f'quantumespresso.{executable.split(".")[0]}',
            prepend_text=prepend_text,
            append_text=append_text
        )
        code.store()
        echo.echo_success(
            f'Code<{code.label}> for {executable} created with pk<{code.pk}> '
            f'on Computer<{computer.label}>'
        )


def _get_code_label(label_template: str, executable: str,
                    computer: orm.Computer) -> Union[Tuple[str, str], Tuple[None, str]]:
    existing_label = None
    label = (
        label_template.replace('{}', executable.replace('.x', '')) if label_template else executable.replace('.x', '')
    )

    # Check if code already exists
    existing_labels = orm.QueryBuilder().append(orm.Computer, filters={
        'label': computer.label
    }, tag='computer').append(
        orm.InstalledCode,
        with_computer='computer',
        filters={
            'label': {
                'like': f'{label}%'
            }
        },
        tag='code',
        project='label'
    ).all(flat=True)
    pattern = re.compile(rf'^{re.escape(label)}(-\d+)?$')
    n_existing_codes = len([l for l in existing_labels if pattern.match(l)])

    if n_existing_codes > 0:
        existing_label = label
        label = f'{label}-{n_existing_codes+1}'

    return existing_label, label


def _get_executable_paths(
    prepend_text: str,
    executables: Union[str, List[str]],
    computer: orm.Computer,
    directory: str,
) -> List[Tuple[str, str]]:
    """Return the absolute paths of the executables on the given computer."""
    user = orm.User.collection.get_default()

    try:
        authinfo = computer.get_authinfo(user)
    except NotExistent:
        echo.echo_critical(f'Computer<{computer.label}> is not yet configured for user<{user.email}>')

    if not authinfo.enabled:
        echo.echo_warning(f'Computer<{computer.label}> is disabled for user<{user.email}>')
        click.confirm('Do you really want to test it?', abort=True)

    transport = authinfo.get_transport()

    executable_paths = []
    with transport:
        for executable in executables:
            if not directory:
                if prepend_text:
                    transport.exec_command_wait(prepend_text)

                # Add prepend text to the command? Sometimes one might need to use modules before finding the executable
                which_ret_val, exec_path, which_stderr = transport.exec_command_wait(f'which {executable}')

                if which_ret_val != 0:
                    echo.echo_error(
                        f'Failed to determine the path of the executable<{executable}> on '
                        f'computer<{computer.label}>:\n\t{which_stderr}'
                    )
                    continue

                exec_path = exec_path.strip()
            else:
                directory = pathlib.PurePosixPath(directory)
                if not directory.is_absolute():
                    echo.echo_error(f'Directory<{directory}> is not an absolute path.')
                    continue
                if not transport.path_exists(directory):
                    echo.echo_error(f'Directory<{directory}> does not exist on computer<{computer.label}>')
                    continue
                exec_path = directory / executable
            # also check that executable is executable and exists
            executable_paths.append((executable, str(exec_path)))

    return executable_paths


# This is just an example of how a test command to test a
# PW code could look like. One should rethink the submission, e.g., which queue
# and resources to be used


@click.argument('codes', nargs=-1, required=True)
def test_code(codes):
    """Test a `code` instance for a Quantum ESPRESSO executable."""
    from aiida.engine import submit
    from aiida.orm import load_code
    from ase.build import bulk

    from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

    # Create a silicon crystal structure
    structure = bulk('Si')
    structure = orm.StructureData(ase=structure)
    assert isinstance(codes, (list, tuple))
    for code_label in codes:
        code = load_code(code_label)
        plugin = code.default_calc_job_plugin.split('.')[-1]
        if plugin != 'pw':
            echo.echo_warning(f'Code<{code.label}> is not a pw.x code, skipping test.')
            continue
        # Create a Quantum ESPRESSO calculation
        builder = PwBaseWorkChain.get_builder_from_protocol(
            code=code,
            structure=structure,
            protocol='fast',
            overrides={'pw': {
                'parameters': {
                    'ELECTRONS': {
                        'electron_maxstep': 2,
                        'scf_must_converge': False
                    }
                }
            }}
        )
        echo.echo_report(f'Submitting PwBaseWorkChain with code<{code.label}>')
        process = submit(builder)
        echo.echo_report(f'PwBaseWorkChain<{process.pk}> submitted')
