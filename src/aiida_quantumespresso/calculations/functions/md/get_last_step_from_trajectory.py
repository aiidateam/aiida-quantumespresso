from aiida import orm
from aiida.engine import calcfunction
import numpy as np

@calcfunction
def get_last_step_from_trajectory(trajectory):

    """The last step of a trajectory is extracted and used as a thermalised trajectory 
    from which a subsequent calculation can be performed. This is useful to separate 
    the equilibration phase from the production phase in an MD simulation.
    """

    if not isinstance(trajectory, orm.TrajectoryData):
        raise Exception('All my inputs have to be instances of TrajectoryData')
    traj = orm.TrajectoryData()
    for arrname in trajectory.get_arraynames():
        if arrname in ('symbols', 'atomic_species_name'):
            traj.set_array(arrname, trajectory.get_array(arrname))
        elif arrname in ('steps', 'times', 'walltimes'):
            traj.set_array(arrname, np.array([trajectory.get_array(arrname)[0]]))
        else:
            traj.set_array(arrname, np.array([trajectory.get_array(arrname)[-1]]))

    [traj.set_attribute(k, v) for k, v in trajectory.attributes_items() if not k.startswith('array|')]
    if 'timestep_in_fs' in trajectory.attributes:
        traj.set_attribute('sim_time_fs', traj.get_array('steps').size * trajectory.get_attribute('timestep_in_fs'))
    return {'last_step_trajectory': traj}