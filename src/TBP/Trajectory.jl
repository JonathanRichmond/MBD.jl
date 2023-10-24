"""
TBP trajectory wrapper

Author: Jonathan Richmond
C: 10/23/23
"""

import MBD: TBPTrajectory

"""
    shallowClone(trajectory)

Return copy of TBP trajectory object

# Arguments
- `trajectory::TBPTrajectory`: TBP trajectory object
"""
function shallowClone(trajectory::TBPTrajectory)
    object = TBPTrajectory(trajectory.initialCondition, trajectory.dynamicsModel)
    object.a = trajectory.a
    object.E = trajectory.E
    object.e = trajectory.e
    object.h = trajectory.h
    object.i = trajectory.i
    object.Omega = trajectory.Omega
    object.omega = trajectory.omega
    object.r_a = trajectory.r_a
    object.r_p = trajectory.r_p
    object.theta = trajectory.theta
    object.TOF = trajectory.TOF

    return object
end
