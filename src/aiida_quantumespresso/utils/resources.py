# -*- coding: utf-8 -*-
"""Utilities for calculation job resources."""
from math import ceil, exp

import numpy as np

from aiida_quantumespresso.utils.defaults.calculation import pw as pw_defaults


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


def get_automatic_parallelization_options(max_num_machines=1, max_wallclock_seconds=1800):  # pylint: disable=invalid-name
    """Return an instance of the automatic parallelization options dictionary.

    :param max_num_machines: set the number of nodes, default=1
    :param max_wallclock_seconds: set the maximum number of wallclock seconds, default=1800
    """
    return {
        'max_num_machines': max_num_machines,
        'target_time_seconds': 0.5 * max_wallclock_seconds,
        'max_wallclock_seconds': max_wallclock_seconds
    }


def get_pw_parallelization_parameters(
    calculation,
    max_num_machines,
    target_time_seconds,
    max_wallclock_seconds,
    calculation_mode='scf',
    round_interval=1800,
    scaling_law=(exp(-16.1951988), 1.22535849)
):
    """Guess optimal choice of parallelzation parameters for a PwCalculation based on a completed initialization run.

    :param calculation: an initial pw calculation (only initialization is sufficient),
        to get number of k-points, of electrons, of spins, fft grids, etc.
    :param max_num_machines: the maximum allowed number of nodes to be used
    :param target_time_seconds: time the calculation should take finally for the user
    :param max_wallclock_seconds: maximum allowed walltime the calculation should take
    :param calculation_mode: kind of calculation_mode to be performed
        ('scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md')
    :param round_interval: the interval in seconds to which the estimated time in the results
        will be rounded up, to determine the max_wallclock_seconds that should be set
    :param scaling_law: list or tuple with 2 numbers giving the
        fit parameters for a power law expressing the single-CPU time to do
        1 scf step, for 1 k-point, 1 spin and 1 small box of the fft grid,
        as a function of number of electrons, in the form: normalized_single_CPU_time = A*n_elec^B
        where A is the first number and B the second.
        Default values were obtained on piz-dora (CSCS) in 2015, on a set of
        4370 calculations (with a very rough fit).

    :return: a dictionary with suggested parallelization parameters with the following keys
        * npools: the number of pools to use in the cmdline setting
        * num_machines: the recommended number of nodes
        * num_mpiprocs_per_machine: the recommended number of processes per nodes
        * estimated_time: the estimated time the calculation should take in seconds
        * max_wallclock_seconds: the recommended max_wall_clock_seconds setting based on the estimated_time value and
        the round_interval argument

    .. note:: If there was an out-of-memory problem during the initial
        calculation, the number of machines is increased.
    """
    # pylint: disable=invalid-name
    from math import gcd

    default_num_mpiprocs_per_machine = calculation.computer.get_default_mpiprocs_per_machine()

    input_parameters = calculation.inputs.parameters.get_dict()
    output_parameters = calculation.outputs.output_parameters.get_dict()
    electron_settings = input_parameters.get('ELECTRONS', {})

    nspin = output_parameters['number_of_spin_components']
    nbands = output_parameters['number_of_bands']
    nkpoints = output_parameters['number_of_k_points']
    nsteps = electron_settings.get('electron_maxstep', pw_defaults.electron_maxstep)
    fft_grid = output_parameters['fft_grid']

    # Determine expected number of scf iterations. In the case of scf-like modes with relax or
    # dynamics steps, we assume an average of 6 steps. All others are single step calculations
    if calculation_mode in ['scf']:
        niterations = nsteps
    elif calculation_mode in ['relax', 'md', 'vc-relax', 'vc-md']:
        niterations = nsteps * 6
    else:
        niterations = 1

    # Compute an estimate single-CPU time
    time_single_cpu = np.prod(fft_grid) * nspin * nkpoints * niterations * scaling_law[0] * nbands**scaling_law[1]

    # The number of nodes is the maximum number we can use that is dividing nkpoints
    num_machines = max([m for m in range(1, max_num_machines + 1) if nkpoints % m == 0])

    # If possible try to make number of kpoints even by changing the number of machines
    if (
        num_machines == 1 and nkpoints > 6 and max_num_machines > 1 and
        time_single_cpu / default_num_mpiprocs_per_machine > target_time_seconds
    ):
        num_machines = max([m for m in range(1, max_num_machines + 1) if (nkpoints + 1) % m == 0])

    # Now we will try to decrease the number of processes per machine (by not more than one fourth)
    # until we manage to get an efficient plane wave parallelization
    # (i.e. number of procs per pool dividing the third dimension of the fft grid)
    num_mpiprocs_per_machine = default_num_mpiprocs_per_machine
    successful = False
    while num_mpiprocs_per_machine >= 0.75 * default_num_mpiprocs_per_machine:
        if nkpoints % num_machines != 0:
            npools = num_machines
        else:
            npools = num_machines * gcd(num_mpiprocs_per_machine, nkpoints / num_machines)
        if fft_grid[2] % num_mpiprocs_per_machine / (npools / num_machines) == 0:
            successful = True
            break
        num_mpiprocs_per_machine -= 1

    if not successful:
        num_mpiprocs_per_machine = default_num_mpiprocs_per_machine
        if nkpoints % num_machines != 0:
            npools = num_machines
        else:
            npools = num_machines * gcd(num_mpiprocs_per_machine, nkpoints / num_machines)

    # Increase the number of machines in case of memory problem during initialization
    if calculation.get_scheduler_stderr() and 'OOM' in calculation.get_scheduler_stderr():
        num_machines = max([i for i in range(num_machines, max_num_machines + 1) if i % num_machines == 0])

    estimated_time = time_single_cpu / (num_mpiprocs_per_machine * num_machines)
    max_wallclock_seconds = min(ceil(estimated_time / round_interval) * round_interval, max_wallclock_seconds)

    result = {
        'resources': {
            'num_machines': num_machines,
            'num_mpiprocs_per_machine': num_mpiprocs_per_machine,
            'tot_num_mpiprocs': num_machines * num_mpiprocs_per_machine,
        },
        'max_wallclock_seconds': max_wallclock_seconds,
        'estimated_time': estimated_time,
        'npools': npools,
    }

    return result
