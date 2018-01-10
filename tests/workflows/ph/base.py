#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm.data.base import Bool, Int
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.utils import CalculationFactory, WorkflowFactory
from aiida.work.run import run
from aiida_quantumespresso.utils.resources import get_default_options

PwCalculation = CalculationFactory('quantumespresso.pw')
PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')

def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the PhBaseWorkChain for a previously completed PwCalculation',
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references QE ph.x'
    )
    parser.add_argument(
        '-p', type=int, required=True, dest='parent_calc',
        help='the node id of the parent PwCalculation'
    )
    parser.add_argument(
        '-q', nargs=3, type=int, default=[2, 2, 2], dest='qpoints', metavar='Q',
        help='define the q-points mesh. (default: %(default)s)'
    )
    parser.add_argument(
        '-m', type=int, default=1, dest='max_num_machines',
        help='the maximum number of machines (nodes) to use for the calculations. (default: %(default)d)'
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
        parent_calc = load_node(args.parent_calc)
    except NotExistent as exception:
        print "Execution failed: failed to load the node for the given parent calculation '{}'".format(args.parent_calc)
        print "Exception report: {}".format(exception)
        return

    if not isinstance(parent_calc, PwCalculation):
        print "The provided parent calculation {} is not of type PwCalculation, aborting...".format(args.parent_calc)
        return

    qpoints = KpointsData()
    qpoints.set_kpoints_mesh(args.qpoints)

    parameters = {
        'INPUTPH': {
        }
    }

    options = get_default_options(args.max_num_machines, args.max_wallclock_seconds)

    inputs = {
        'code': code,
        'qpoints': qpoints,
        'parent_calc': parent_calc,
        'parameters': ParameterData(dict=parameters),
        'options': ParameterData(dict=options),
    }

    run(PhBaseWorkChain, **inputs)


def main():
    """
    Setup the parser to retrieve the command line arguments and pass them to the main execution function.
    """
    parser = parser_setup()
    args = parser.parse_args()
    result = execute(args)


if __name__ == "__main__":
    main()