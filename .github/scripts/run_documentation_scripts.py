#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run scripts from the documentation to verify they run without exceptions."""
import pathlib
import subprocess

from aiida.cmdline.utils.decorators import with_dbenv
import click

FILEPATH_DOCS = pathlib.Path(__file__).parent.parent.parent / 'docs' / 'source'

# Define the scripts that should be run and in which order
FILEPATH_SCRIPTS = (
    FILEPATH_DOCS / 'tutorials' / 'include' / 'scripts' / 'run_pw_basic.py',
    FILEPATH_DOCS / 'tutorials' / 'include' / 'scripts' / 'run_ph_basic.py',
)


def get_pw_calculation_pk():
    """Return the pk of the first ``PwCalculation`` that is found."""
    from aiida.orm import CalcJobNode, QueryBuilder
    return QueryBuilder().append(
        CalcJobNode, filters={
            'process_type': 'aiida.calculations:quantumespresso.pw'
        }, project='id'
    ).first(flat=True)


@with_dbenv()
def main():
    """Run scripts from the documentation to verify they run without exceptions."""

    # Mapping of placeholder strings in the scripts that should be replaced with the return value of the callable
    placeholder_replacement = {
        'IDENTIFIER_PW_CALCULATION': get_pw_calculation_pk,
    }

    for filepath in FILEPATH_SCRIPTS:

        filepath_content = filepath.read_text()

        for key, function in placeholder_replacement.items():
            if key in filepath_content:
                filepath_content.replace(key, str(function()))

        filepath.write_text(filepath_content)

        click.echo(f'Running: {filepath.relative_to(FILEPATH_DOCS)} ...', nl=False)
        try:
            stdout = subprocess.check_output(['verdi', 'run', str(filepath)])
        except subprocess.CalledProcessError:
            click.secho('[ERROR]', fg='red', bold=True)
        else:
            click.secho('[OK]', fg='green', bold=True)
            click.echo(stdout)


if __name__ == '__main__':
    main()
