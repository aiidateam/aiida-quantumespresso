from aiida import orm
from aiida.engine import calcfunction
import numpy as np

SINGULAR_TRAJ_KEYS = ('symbols', 'atomic_species_name')

@calcfunction
def concatenate_trajectory(**kwargs):
    remove_repeated_last_step = kwargs.pop('remove_repeated_last_step')
    for k, v in kwargs.items():
        if not isinstance(v, orm.TrajectoryData):
            raise Exception('All my inputs have to be instances of TrajectoryData')
    sorted_trajectories = list(zip(*sorted(kwargs.items())))[1]
    # Assumming they store the same arrays
    arraynames = sorted_trajectories[0].get_arraynames()
    traj = orm.TrajectoryData()
    for arrname in arraynames:
        if arrname in SINGULAR_TRAJ_KEYS:
            traj.set_array(arrname, sorted_trajectories[0].get_array(arrname))
        else:
            # concatenate arrays
            if len(sorted_trajectories) > 1:
                if remove_repeated_last_step:
                    # remove last step that is repeated when restarting, keep the very last
                    # typically used for pinball binary
                    traj.set_array(arrname, np.concatenate([
                        np.concatenate([t.get_array(arrname)[:-1] for t in sorted_trajectories[:-1]]),
                        sorted_trajectories[-1].get_array(arrname)]))
                else:
                    # just concatenate, typically used for vanilla QE
                    traj.set_array(arrname, np.concatenate([t.get_array(arrname) for t in sorted_trajectories]))
            else:
                traj.set_array(arrname, np.concatenate([t.get_array(arrname) for t in sorted_trajectories]))
    [traj.set_attribute(k, v) for k, v in sorted_trajectories[0].attributes_items() if not k.startswith('array|')]
    # TODO: if the timestep is different, then stride the denser trajectory so that each step is at the same 
    # level, and they also both must have the same strides
    if 'timestep_in_fs' in sorted_trajectories[0].attributes:
        traj.set_attribute('sim_time_fs', traj.get_array('steps').size * sorted_trajectories[0].get_attribute('timestep_in_fs'))
    return {'concatenated_trajectory': traj}