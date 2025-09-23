# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Test the setup tools."""
from contextlib import nullcontext

import pytest

from aiida_quantumespresso.tools.code_setup import get_executable_paths


@pytest.fixture
def create_pw_executable(tmp_path):
    """Create a dummy executable `pw.x` in a temp dir and return its path."""

    def _factory(subdir='bin'):
        executable_path = tmp_path / subdir / 'pw.x'
        executable_path.parent.mkdir(parents=True, exist_ok=True)
        executable_path.touch()
        executable_path.chmod(0o755)  # Needs to be executable for `which` to find it
        return executable_path

    return _factory


@pytest.mark.parametrize(
    'prepend_text,error,error_message', (
        ('export PATH={path}:$PATH', None, ''),
        ('echo lala\nexport PATH={path}:$PATH', None, ''),
        ('# Comment\nexport PATH={path}:$PATH', None, ''),
        ('', FileNotFoundError, 'Error: the `which` command returned an empty output.'),
        ('echo lala', FileNotFoundError, 'Error: the `which` command returned an empty output.'),
        ('expoat PATH={path}:$PATH', FileNotFoundError, 'expoat: command not found'),
        ('export PATH={path}WRONG:$PATH', FileNotFoundError, 'Error: the `which` command returned an empty output.'),
    )
)
def test_get_executable_paths_prepend_text(create_pw_executable, fixture_localhost, prepend_text, error, error_message):
    """Tests the `get_executable_paths` function for various prepend texts."""
    pw_executable = create_pw_executable()
    prepend_text = prepend_text.format(path=pw_executable.parent.as_posix())

    # `pytest.raises` does not allow the error to be `None`, so we need to adapt the context
    context = pytest.raises(error, match=error_message) if error else nullcontext()

    with context:
        result = get_executable_paths(
            executable_tuple=(pw_executable.name,),
            computer=fixture_localhost,
            prepend_text=prepend_text,
        )
        assert result == {pw_executable.name: pw_executable.as_posix()}


def test_get_executable_paths_with_quoted_path(create_pw_executable, fixture_localhost):
    """Tests the `get_executable_paths` function for quotes path that has a space."""
    pw_executable = create_pw_executable('this-dir-has-a space/bin')
    prepend_text = f'export PATH="{pw_executable.parent.as_posix()}":$PATH'

    result = get_executable_paths(
        executable_tuple=(pw_executable.name,),
        computer=fixture_localhost,
        prepend_text=prepend_text,
    )
    assert result == {pw_executable.name: pw_executable.as_posix()}
