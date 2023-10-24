"""
TBP trajectory wrapper

Author: Jonathan Richmond
C: 10/23/23
"""

import MBD: TBPTrajectory

export getCartesianState

"""
    getCartesianState(trajectory, theta)

Return 2BP trajectory object with Cartesian state

# Arguments
- `trajectory::TBPTrajectory`: TBP trajectory object
- `theta::Float64`: True anomaly [rad]
"""
function getCartesianState(trajectory::TBPTrajectory, theta::Float64)
    eccentricAnomaly::Float64 = atan(sqrt(1-trajectory.e^2)*sin(theta), trajectory.e+cos(theta))
    circularRadius::Float64 = trajectory.a*(1-trajectory.e*cos(eccentricAnomaly))
    r_E::Vector{Float64} = [circularRadius*cos(theta), circularRadius*sin(theta), 0]
    v_E::Vector{Float64} = (sqrt(trajectory.dynamicsModel.systemData.gravParam*trajectory.a)/circularRadius).*[-sin(eccentricAnomaly), sqrt(1-trajectory.e^2)*cos(eccentricAnomaly), 0]
    C::Matrix{Float64} = [cos(trajectory.omega)*cos(trajectory.Omega)-sin(trajectory.omega)*cos(trajectory.i)*sin(trajectory.Omega) -sin(trajectory.omega)*cos(trajectory.Omega)-cos(trajectory.omega)*cos(trajectory.i)*sin(trajectory.Omega) 0; cos(trajectory.omega)*sin(trajectory.Omega)+sin(trajectory.omega)*cos(trajectory.i)*cos(trajectory.Omega) -sin(trajectory.omega)*sin(trajectory.Omega)+cos(trajectory.omega)*cos(trajectory.i)*cos(trajectory.Omega) 0; sin(trajectory.omega)*sin(trajectory.i) cos(trajectory.omega)*sin(trajectory.i) 0]
    pos_dim::Vector{Float64} = C*r_E
    vel_dim::Vector{Float64} = C*v_E
    stateTrajectory = TBPTrajectory(append!(pos_dim, vel_dim), trajectory.dynamicsModel)
    angularMomentum::Vector{Float64} = LinearAlgebra.cross(pos_dim, vel_dim)
    stateTrajectory.h = LinearAlgebra.norm(angularMomentum)
    stateTrajectory.i = trajectory.i
    stateTrajectory.Omega = trajector.Omega
    stateTrajectory.e = trajectory.e
    stateTrajectory.a = trajectory.a
    stateTrajectory.r_p = trajectory.a*(1-trajectory.e)
    stateTrajectory.r_a = trajectory.a*(1+trajectory.e)
    stateTrajectory.omega = trajectory.omega
    stateTrajectory.theta = theta
    stateTrajectory.E = -trajectory.dynamicsModel.systemData.gravParam/(2*trajectory.a)

    return stateTrajectory
end

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
