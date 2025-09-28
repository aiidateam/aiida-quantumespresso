# -*- coding: utf-8 -*-
"""Module for setting up AiiDA components for Quantum ESPRESSO."""

from pathlib import PurePosixPath
import re
from typing import Union
import warnings

from aiida import orm


def setup_codes(
    computer: orm.Computer,
    executables: Union[str, list[str], tuple[str]],
    directory: Union[str, None] = None,
    label_template: str = '',
    prepend_text: str = '',
    append_text: str = '',
    increment_label: bool = False,
):
    """Automatically create `orm.Code` instances for Quantum ESPRESSO executables.

    Specify the target `computer` and a single executable or a list of Quantum ESPRESSO `executables` to create codes
    for. You can provide multiple executables as a list or tuple.

    Use `directory` to point to the directory containing the executables. If not provided, the command will try to find
    the executables in the `PATH` on the (remote) computer. Consider adding a `prepend_text` if modules need to be
    loaded to find the executables.

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
    :param increment_label: if a code with the same label already exists on the computer, append an incremental suffix
        to the label, e.g. `pw-2` if `pw` already exists. Else a warning will be raised but no code is created.
    """
    executable_tuple = (executables,) if isinstance(executables, str) else tuple(executables)

    created = []

    for executable, exec_path in get_executable_paths(executable_tuple, computer, prepend_text, directory).items():
        existing_label, label = get_code_label(label_template=label_template, executable=executable, computer=computer)

        if existing_label:
            if increment_label:
                warnings.warn(
                    f'Code with label "{existing_label}" already exists on computer "{computer.label}". '
                    f'`increment_label` is True: creating another one with label "{label}".'
                )
            else:
                warnings.warn(
                    f'Code with label "{existing_label}" already exists on computer "{computer.label}". '
                    f'Skipping executable "{executable}". Consider setting `increment_label=True` or '
                    'specifying a new `label_template`.'
                )
                continue

        code = orm.InstalledCode(
            label=label,
            computer=computer,
            filepath_executable=exec_path,
            default_calc_job_plugin=f"quantumespresso.{executable.split('.')[0]}",
            prepend_text=prepend_text,
            append_text=append_text,
        )
        code.store()
        created.append(code)

    return created


def get_executable_paths(
    executable_tuple: tuple[str],
    computer: orm.Computer,
    prepend_text: str = '',
    directory: Union[str, None] = None,
) -> dict:
    """Return a mapping from executable names to absolute paths on the given computer.

    If `directory` is provided, each path is constructed as `directory`/`executable`.
    Otherwise, the `prepend_text` is combined with `which` to locate executables in the PATH.
    """
    executable_paths = {}

    with computer.get_transport() as transport:
        for executable in executable_tuple:
            if directory is None:
                combined_prepend_text = '\n'.join((computer.get_prepend_text(), prepend_text))
                return_value, stdout, stderr = transport.exec_command_wait(
                    command=f'. /dev/stdin > /dev/null && which {executable}', stdin=combined_prepend_text
                )

                if return_value != 0 or not stdout.strip():
                    msg = f'Failed to determine the path of executable<{executable}> on computer<{computer.label}>.\n'
                    if stderr:
                        msg += f'Error: {stderr}'
                    elif not stdout.strip():
                        msg += 'Error: the `which` command returned an empty output.\n'
                    msg += (
                        '\nDouble-check the `prepend_text` and executables and/or specify the full path with the '
                        '`directory` input.'
                    )
                    raise FileNotFoundError(msg)

                executable_paths[executable] = PurePosixPath(stdout.strip()).as_posix()
            else:
                directory = PurePosixPath(directory)

                if not directory.is_absolute():
                    raise ValueError(f'Directory<{directory}> is not an absolute path.')

                exec_path = (directory / executable).as_posix()

                if not transport.path_exists(exec_path):
                    warnings.warn(f'Could not find executable<{exec_path}> on computer<{computer.label}>')
                executable_paths[executable] = exec_path

    return executable_paths


def get_code_label(
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
