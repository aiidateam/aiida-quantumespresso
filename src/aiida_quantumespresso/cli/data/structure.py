# -*- coding: utf-8 -*-
"""Command line utilities to create and inspect `StructureData` nodes."""
from aiida.cmdline.params import options
from aiida.cmdline.utils import decorators, echo
import click

from . import cmd_data


@cmd_data.group('structure')
def cmd_structure():
    """Commands to create and inspect `StructureData` nodes."""


@cmd_structure.command('import')
@click.argument('filename', type=click.File('r'))
@options.DRY_RUN()
@decorators.with_dbenv()
def cmd_import(filename, dry_run):
    """Import a `StructureData` from a Quantum ESPRESSO input file."""
    from aiida_quantumespresso.tools.pwinputparser import PwInputFile

    with open(filename, 'r', encoding='utf-8') as input_file:
        parser = PwInputFile(input_file.read())
    structure = parser.get_structuredata()
    formula = structure.get_formula()

    if dry_run:
        echo.echo_success(f'parsed structure with formula {formula}')
    else:
        structure.store()
        echo.echo_success(f'parsed and stored StructureData<{structure.pk}> with formula {formula}')
