# -*- coding: utf-8 -*-
"""Utilities for calculation job resources."""


def create_scheduler_resources(scheduler, base, goal):
    """Return a dictionary with scheduler resources.

    This function will take a dictionary 'base', which contains scheduler resource settings that can be set for the
    'resources' key in the 'options' input of a CalcJob, and update it with target scheduler resource value from the
    dictionary 'goal', checking that the resource settings are allowed for the specific scheduler implementation. The
    function should than return a dictionary that can be used directly in the 'options' input of the CalcJob. Any keys
    in either 'base' or 'goal' that are not recognized as valid resource keys by the scheduler class, will be removed
    and will be kept otherwise. For settings that appear both in 'base' and 'goal', the value in 'goal' will be used.

    :param scheduler: instance of the Scheduler class
    :param base: dictionary with base scheduler resource settings
    :param goal: dictionary with target scheduler resource settings
    :returns: dictionary that can be used for the 'resources' key in the 'options' dictionary of a CalcJob
    """
    base.update(goal)

    resources = {}
    for key, value in base.items():
        if key in scheduler._job_resource_class.get_valid_keys():  # pylint: disable=protected-access
            resources[key] = value

    try:
        job_resource = scheduler.create_job_resource(**resources)
    except TypeError as exception:
        raise ValueError(f'failed to create job resources for {scheduler.__class__} scheduler') from exception

    return {key: value for key, value in job_resource.items() if value is not None}


def cmdline_remove_npools(cmdline):
    """Remove all options related to npools in the `settings.cmdline` input.

    The cmdline setting is a list of strings that typically looks something like:

        cmdline = ['-nk', '4', '-ntg', '8']

    This function will remove all occurrences of '-nk', '-npool', '-npools', which are
    all synonymous flags, and the directly following element, which should be the integer

    :param cmdline: the cmdline setting which is a list of string directives
    :return: the new cmdline setting
    """
    return [
        e for i, e in enumerate(cmdline)
        if (e not in ('-npools', '-npool', '-nk') and cmdline[i - 1] not in ('-npools', '-npool', '-nk'))
    ]


def get_default_options(max_num_machines=1, max_wallclock_seconds=1800, with_mpi=False):
    """Return an instance of the options dictionary with the minimally required parameters for a `CalcJob`.

    :param max_num_machines: set the number of nodes, default=1
    :param max_wallclock_seconds: set the maximum number of wallclock seconds, default=1800
    :param with_mpi: whether to run the calculation with MPI enabled
    """
    return {
        'resources': {
            'num_machines': int(max_num_machines)
        },
        'max_wallclock_seconds': int(max_wallclock_seconds),
        'withmpi': with_mpi,
    }
