# -*- coding: utf-8 -*-
"""Command line utilities to create and inspect `StructureData` nodes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options
from aiida.cmdline.utils import decorators, echo

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

    parser = PwInputFile(filename)
    structure = parser.get_structuredata()
    formula = structure.get_formula()

    if dry_run:
        echo.echo_success('parsed structure with formula {}'.format(formula))
    else:
        structure.store()
        echo.echo_success('parsed and stored StructureData<{}> with formula {}'.format(structure.pk, formula))
