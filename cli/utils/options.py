# -*- coding: utf-8 -*-
from __future__ import absolute_import 
import click
from cli.utils import validators

class overridable_option(object):
    """
    Wrapper around click option that allows to store the name
    and some defaults but also to override them later, for example
    to change the help message for a certain command.
    """

    def __init__(self, *args, **kwargs):
        """
        Store the defaults.
        """
        self.args = args
        self.kwargs = kwargs

    def __call__(self, **kwargs):
        """
        Override kwargs (no name changes) and return option
        """
        kw_copy = self.kwargs.copy()
        kw_copy.update(kwargs)
        return click.option(*self.args, **kw_copy)

code = overridable_option(
    '-c', '--code', type=click.STRING, required=True,
    callback=validators.validate_code,
    help='the label of the AiiDA code object to use'
)

structure = overridable_option(
    '-s', '--structure', type=click.INT, required=True,
    callback=validators.validate_structure,
    help='the node id of the structure'
)

pseudo_family = overridable_option(
    '-p', '--pseudo-family', type=click.STRING, required=True,
    callback=validators.validate_pseudo_family,
    help='the name of the pseudo potential family to use'
)

kpoint_mesh = overridable_option(
    '-k', '--kpoint-mesh', 'kpoints', nargs=3, type=click.INT, default=[2, 2, 2], show_default=True,
    callback=validators.validate_kpoint_mesh,
    help='the number of points in the kpoint mesh along each basis vector'
)

max_num_machines = overridable_option(
    '-m', '--max-num-machines', type=click.INT, default=1, show_default=True,
    help='the maximum number of machines (nodes) to use for the calculations'
)

max_wallclock_seconds = overridable_option(
    '-w', '--max-wallclock-seconds', type=click.INT, default=1800, show_default=True,
    help='the maximum wallclock time in seconds to set for the calculations'
)

automatic_parallelization = overridable_option(
    '-a', '--automatic-parallelization', is_flag=True, default=False, show_default=True,
    help='enable the automatic parallelization option of the workchain'
)

clean_workdir = overridable_option(
    '-x', '--clean-workdir', is_flag=True, default=False, show_default=True,
    help='clean the remote folder of all the launched calculations after completion of the workchain'
)