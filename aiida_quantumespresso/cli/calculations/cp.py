"""
Check that the CP calculation might not work.
"""

from __future__ import absolute_import, print_function

from pprint import pprint

import click
from aiida.cmdline.params import options, types
from aiida.engine import run
from aiida.orm import Dict, StructureData, QueryBuilder
from aiida.orm.nodes.data.upf import get_pseudos_from_structure
from ase.spacegroup import crystal

from aiida_quantumespresso.cli.utils import options as options_qe


def silicon_structure():
    """
    Builds a new silicon structure.
    """
    label = '__example__{name}'.format(name=__name__)
    structure = QueryBuilder().append(StructureData, filters={'label': label}).first()
    if not structure:
        alat = 5.4
        ase_structure = crystal(
            "Si",
            [(0, 0, 0)],
            spacegroup=227,
            cellpar=[alat, alat, alat, 90, 90, 90],
            primitive_cell=True,
        )
        structure = StructureData(ase=ase_structure)
        structure.label = label
        structure.store()
    else:
        # .frist() returns none if it doesn't exist or a list of one item if objects exist
        structure = structure[0]
    return structure.uuid


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.cp'))
@options_qe.STRUCTURE(default=silicon_structure)
@options_qe.PSEUDO_FAMILY(required=True)
@options.DRY_RUN()
def cli(code, structure, pseudo_family, dry_run):
    """
    Little check on the cp calculation thing.
    """

    builder = code.get_builder()

    builder.structure = structure

    # Parameters that we of course don't know about
    parameters = Dict(
        dict={
            'CONTROL': {
                'calculation': "cp",
                'restart_mode': "from_scratch",
                'wf_collect': False,
                'iprint': 1,
                'isave': 100,
                'dt': 3.0,
                'max_seconds': 25 * 60,
                'nstep': 10,
            },
            'SYSTEM': {
                'ecutwfc': 30.0,
                'ecutrho': 240.0,
                'nr1b': 24,
                'nr2b': 24,
                'nr3b': 24,
            },
            'ELECTRONS': {
                'electron_damping': 1.0e-1,
                'electron_dynamics': "damp",
                'emass': 400.0,
                'emass_cutoff': 3.0,
            },
            'IONS': {
                'ion_dynamics': "none"
            },
        })
    builder.parameters = parameters

    pseudos = get_pseudos_from_structure(structure, pseudo_family)
    builder.pseudos = pseudos

    builder.metadata.dry_run = dry_run
    builder.metadata.store_provenance = not dry_run
    builder.metadata.options = {"resources": {"num_machines": 1}}
    results, node = run.get_node(builder)
    print('The calculation node is:')
    pprint(node)

    print('\nAnd the results are:')
    pprint(results)
