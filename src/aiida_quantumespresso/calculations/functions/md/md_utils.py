# Utility functions for MD workchains that are not supposed to be calcfunctions

from aiida import orm
from aiida.common.links import LinkType


def get_completed_number_of_steps(calc):
    """Read the number of steps from the trajectory."""
    traj = calc.outputs.output_trajectory
    nstep = (
        calc.inputs.parameters.get_attribute('CONTROL').get('iprint', 1) * (traj.get_attribute('array|positions')[0])
    )
    return nstep


def get_total_trajectory(workchain, previous_trajectory=None, store=False, remove_repeated_last_step=False):
    """Collect all the trajectory segment and concatenate them."""

    from aiida_quantumespresso.calculations.functions.md.concatenate_trajectory import concatenate_trajectory

    qb = orm.QueryBuilder()
    qb.append(orm.WorkChainNode, filters={'uuid': workchain.uuid}, tag='replay')
    # TODO: Are filters on the state of the calculation needed here?
    # TODO: add project on extras.discard_trajectory, traj_d defined to skip them
    qb.append(
        orm.CalcJobNode,
        with_incoming='replay',
        edge_filters={'type': LinkType.CALL_CALC.value, 'label': {'like': 'iteration_%'}},
        edge_project='label',
        tag='calc',
        edge_tag='rc',
    )
    qb.append(
        orm.TrajectoryData, with_incoming='calc', edge_filters={'label': 'output_trajectory'}, project=['*'], tag='traj'
    )
    traj_d = {
        item['rc']['label'].replace('iteration_', 'trajectory_'): item['traj']['*'] for item in qb.iterdict()
    }  ## if not extras.discard_trajectory

    # adding the trajectory of previous MD run, if it exists
    if previous_trajectory:
        traj_d.update({'trajectory_00': previous_trajectory})

    # if I have produced several trajectories, I concatenate them here: (no need to sort them)
    if len(traj_d) > 1:
        traj_d['metadata'] = {'call_link_label': 'concatenate_trajectory', 'store_provenance': store}
        # If I want to start the next MD from the position and velocities of the previous one,
        # i.e. using `previous_trajectory`, I need to ensure that I do not duplicate
        # the last step from the previous trajectory.
        # In case of restarting from wavefunctions this must be false
        traj_d.update({'remove_repeated_last_step': orm.Bool(remove_repeated_last_step)})
        res = concatenate_trajectory(**traj_d)
        return res['concatenated_trajectory']
    elif len(traj_d) == 1:
        # no reason to concatenate if I have only one trajectory (saves space in repository)
        return list(traj_d.values())[0]
    else:
        return None
