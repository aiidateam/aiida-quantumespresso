#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm.data.base import Str
from aiida.orm.data.upf import UpfData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.utils import WorkflowFactory
from aiida.work.run import run

PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')

def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description="""
        Run the PwBandStructureWorkChain for a given input structure 
        to compute the band structure for the relaxed structure
        """,
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
        '-m', type=str, dest='protocol', default='standard',
        help='the protocol to use for the calculation  (default: %(default)s)'
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
        pseudo_family = UpfData.get_upf_group(args.pseudo_family)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the pseudo family group '{}'".format(args.pseudo_family)
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

    run(
        PwBandStructureWorkChain,
        code=code,
        structure=structure,
        pseudo_family=Str(args.pseudo_family),
        protocol=Str(args.protocol)
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