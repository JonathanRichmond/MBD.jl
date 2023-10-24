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
    eccentricAnomaly::Float64 = atan(sqrt(1-trajectory.e^2)*sin(theta), e+cos(theta))
    circularRadius::Float64 = trajectory.a*(1-trajectory.e*cos(eccentricAnomaly))
    r0::Vector{Float64} = [circularRadius*cos(theta), circularRadius*sin(theta), 0]
    v0::Vector{Float64} = (sqrt(dynamicsModel.systemData.gravParam*trajectory.a)/circularRadius).*[-sin(eccentricAnomaly), sqrt(1-trajectory.e^2)*cos(eccentricAnomaly), 0]
    C::Matrix{Float64} = [cos(trajectory.omega)*cos(trajectory.Omega)-sin(trajectory.omega)*cos(trajectory.i)*sin(trajectory.Omega) -sin(trajectory.omega)*cos(trajectory.Omega)-cos(trajectory.omega)*cos(trajectory.i)*sin(trajectory.Omega) 0; cos(trajectory.omega)*sin(trajectory.Omega)+sin(trajectory.omega)*cos(trajectory.i)*cos(trajectory.Omega) -sin(trajectory.omega)*sin(trajectory.Omega)+cos(trajectory.omega)*cos(trajectory.i)*cos(trajectory.Omega) 0; sin(trajectory.omega)*sin(trajectory.i) cos(trajectory.omega)*sin(trajectory.i) 0]
    r::Vector{Float64} = C*r0
    v::Vector{Float64} = C*v0

    return getOsculatingOrbitalElements(trajectory.dynamicsModel, append!(r, v))
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
