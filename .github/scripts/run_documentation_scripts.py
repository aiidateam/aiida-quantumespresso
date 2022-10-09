#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run scripts from the documentation to verify they run without exceptions."""
import pathlib
import subprocess

import click

EXCLUDE_FILES = ('conf.py')


def main():
    """Run scripts from the documentation to verify they run without exceptions."""
    filepath_docs = pathlib.Path(__file__).parent.parent.parent / 'docs' / 'source'
    filepath_scripts = filepath_docs.glob('**/*.py')

    for filepath in filter(lambda p: str(p.relative_to(filepath_docs)) not in EXCLUDE_FILES, filepath_scripts):
        click.echo(f'Running: {filepath.relative_to(filepath_docs)} ...', nl=False)

        try:
            stdout = subprocess.check_output(['verdi', 'run', str(filepath)])
        except subprocess.CalledProcessError:
            click.secho('[ERROR]', fg='red', bold=True)
            click.echo(stdout)
        else:
            click.secho('[OK]', fg='green', bold=True)


if __name__ == '__main__':
    main()
