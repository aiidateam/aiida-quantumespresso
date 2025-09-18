# -*- coding: utf-8 -*-
"""Module for setting up `Code` instances for Quantum ESPRESSO executables."""

import pathlib
import re
from typing import Union

from aiida import orm
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent
import click


def create_codes(
    computer: orm.Computer,
    executables: Union[str, list[str], tuple[str]],
    directory: Union[str, None] = None,
    label_template: str = '',
    prepend_text: Union[str, None] = None,
    append_text: Union[str, None] = None,
    on_conflict: str = 'increment',
):
    """Automatically create `orm.Code` instances for Quantum ESPRESSO executables.

    :param computer: the `orm.Computer` on which the executables are located
    :param executables: Quantum ESPRESSO executable(s) to create code(s) for, e.g. `pw.x` or `['pw.x', 'projwfc.x']`
    :param directory: absolute path to a directory where the executable(s) is (are) located
                      If not provided, the command will try to find the executables in the `PATH`
                      on the (remote) computer
    :param label_template: label template for the code instance. Use curly brackets to reference
                           the executable label, e.g. `qe-{}` will create a `Code` with label
                           `qe-pw` for `pw.x`. Defaults to the executable name without `.x` suffix.,
                           e.g. `pw` for `pw.x`.
    :param prepend_text: text to prepend to the execution command, e.g. to load modules
    :param append_text: text to append to the execution command
    :param on_conflict: action to take if a code with the same label already exists on the computer. Can be one of:
                        - 'increment': append an incremental suffix to the label, e.g. `pw-2` if `pw` already exists
                        - 'skip': skip creating the code for that executable
                        - 'prompt': prompt the user to confirm whether to create another code with an incremented
                          label or skip creating the code for that executable
    """
    if isinstance(executables, str):
        executables = [executables]

    return _create_codes(
        computer=computer,
        executables=executables,
        directory=directory,
        label_template=label_template,
        prepend_text=prepend_text,
        append_text=append_text,
        on_conflict=on_conflict,
    )


def _create_codes(
    *,
    computer: orm.Computer,
    executables: list[str],
    directory: Union[str, None] = None,
    label_template: str = '',
    prepend_text: Union[str, None] = None,
    append_text: Union[str, None] = None,
    on_conflict: str = 'increment',
):
    """Automatically create `orm.Code` instances for Quantum ESPRESSO executables."""
    if on_conflict not in {'increment', 'skip', 'prompt'}:
        raise ValueError("on_conflict must be one of 'increment', 'skip', or 'prompt'")

    exec_paths = _get_executable_paths(
        prepend_text=prepend_text, executables=executables, computer=computer, directory=directory
    )

    created = []
    for executable, exec_path in exec_paths:
        existing_label, label = _get_code_label(label_template=label_template, executable=executable, computer=computer)
        if existing_label:
            if on_conflict == 'prompt':
                echo.echo_warning(f'Code with label<{existing_label}> already exists on Computer<{computer.label}>.')
                if not click.confirm(f'Do you want to add another instance with label {label}?'):
                    continue
            elif on_conflict == 'increment':
                echo.echo_info(
                    f'Code with label<{existing_label}> already exists on Computer<{computer.label}>; '
                    f'creating another one with label<{label}>.'
                )
            else:
                echo.echo_warning(f'Skipping executable<{executable}>.')
                echo.echo_report(
                    f'Code with label<{existing_label}> already exists on Computer<{computer.label}>.'
                    "\nSet `on_conflict='increment'` to automatically increment the label to "
                    f'<{label}> or manually change the `label_template`.\n'
                )
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
            f'Code<{code.label}> for {executable} created with pk<{code.pk}> '
            f'on Computer<{computer.label}>'
        )
        created.append(code)

    return created


def _get_executable_paths(
    prepend_text: str,
    executables: Union[str, list[str]],
    computer: orm.Computer,
    directory: str,
) -> list[tuple[str, str]]:
    """Return the absolute paths of the executables on the given computer.

    If `directory` is provided, the path is constructed as `directory`/`executable`.
    If not, the `which` command is used to find the executable in the `PATH`.
    """
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

            executable_paths.append((executable, str(exec_path)))

    return executable_paths


def _get_code_label(
    label_template: str,
    executable: str,
    computer: orm.Computer,
) -> Union[tuple[str, str], tuple[None, str]]:
    """Return a label for the code based on the executable and computer.

    If ``increment_if_exists`` is True and a code with the same base label exists on the computer,
    append an incremental suffix ``-N`` to the label and return the original conflicting label as ``existing_label``.
    If False, return the base label and set ``existing_label`` when a conflict is detected without incrementing.
    """
    existing_label = None
    if label_template:
        label = label_template.replace('{}', executable.replace('.x', ''))
    else:
        label = executable.replace('.x', '')

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
    n_existing_codes = len([l for l in list(existing_labels) if pattern.match(l)])

    if n_existing_codes > 0:
        existing_label = label
        label = f'{label}-{n_existing_codes+1}'

    return existing_label, label
