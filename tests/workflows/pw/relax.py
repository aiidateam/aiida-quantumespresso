#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm.data.base import Str
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.utils import WorkflowFactory
from aiida.work.run import run

PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')

def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the PwRelaxWorkChain for a given input structure',
    )
    parser.add_argument(
        '-k', nargs=3, type=int, default=[2, 2, 2], dest='kpoints', metavar='Q',
        help='define the q-points mesh. (default: %(default)s)'
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references QE pw.x'
    )
    parser.add_argument(
        '-p', type=str, required=True, dest='pseudo_family',
        help='the name of pseudo family to use'
    )
    parser.add_argument(
        '-s', type=int, required=True, dest='structure',
        help='the node id of the structure'
    )
    parser.add_argument(
        '-w', type=int, default=1800, dest='max_wallclock_seconds',
        help='the maximum wallclock time in seconds to set for the calculations. (default: %(default)d)'
    )

    return parser


def execute(args):
    """
    The main execution of the script, which will run some preliminary checks on the command
    line arguments before passing them to the workchain and running it
    """
    try:
        code = Code.get_from_string(args.codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.codename)
        print "Exception report: {}".format(exception)
        return

    try:
        structure = load_node(args.structure)
    except NotExistent as exception:
        print "Execution failed: failed to load the node for the given structure pk '{}'".format(args.structure)
        print "Exception report: {}".format(exception)
        return

    if not isinstance(structure, StructureData):
        print "The provided pk {} for the structure does not correspond to StructureData, aborting...".format(args.parent_calc)
        return

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(args.kpoints)

    parameters = {
        'CONTROL': {
            'restart_mode': 'from_scratch',
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'ecutrho': 240.,
        },
    }
    settings = {}
    options  = {
        'resources': {
            'num_machines': 1,
            'tot_num_mpiprocs': 8,
        },
        'max_wallclock_seconds': args.max_wallclock_seconds,
    }

    run(
        PwRelaxWorkChain,
        code=code,
        structure=structure,
        pseudo_family=Str(args.pseudo_family),
        kpoints=kpoints,
        parameters=ParameterData(dict=parameters),
        settings=ParameterData(dict=settings),
        options=ParameterData(dict=options)
    )


def main():
    """
    Setup the parser to retrieve the command line arguments and pass them to the main execution function.
    """
    parser = parser_setup()
    args   = parser.parse_args()
    result = execute(args)


if __name__ == "__main__":
    main()