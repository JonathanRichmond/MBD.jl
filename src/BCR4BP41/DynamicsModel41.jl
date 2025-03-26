"""
BCR4BP P4-B1 dynamics model wrapper

Author: Jonathan Richmond
C: 2/20/25
U: 3/26/25
"""

import StaticArrays
import MBD: BCR4BP41DynamicsModel

export appendExtraInitialConditions, checkSTM, evaluateEquations, getEpochDependencies
export getEpochTime, getEquationsOfMotion, getExcursion, getHamiltonian, getParameterDependencies
export getPrimaryState, getPseudopotentialJacobian, getStateSize, getStateTransitionMatrix
export gettheta2, get12MassRatio, get2BApproximation, get4Distance, get4Mass, get41CharLength
export get41CharTime, get41MassRatio, isEpochIndependent, primaryEclipticToRotating41
export rotating41ToPrimaryEcliptic, rotating41ToRotating12

"""
    appendExtraInitialConditions(dynamicsModel, q0_simple, outputEquationType)

Return state vector with extra initial conditions

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial state vector [ndim]
- `outputEquationType::EquationType`: Output state EOM type
"""
function appendExtraInitialConditions(dynamicsModel::BCR4BP41DynamicsModel, q0_simple::Vector{Float64}, outputEquationType::MBD.EquationType)
    n_in::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    (Int16(length(q0_simple)) == n_in) || throw(ArgumentError("State vector length is $(length(q0_simple)), but should be $n_in"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    n_STM::Int16 = getStateSize(dynamicsModel, MBD.STM)
    n_out::Int16 = getStateSize(dynamicsModel, outputEquationType)
    q0_out::Vector{Float64} = zeros(Float64, n_out)
    if n_in >= n_out
        q0_out = q0_simple[1:n_out]
    else
        q0_out[1:n_in] = q0_simple
        [q0_out[i] = 1 for i in n_simple+1:n_simple+1:n_STM]
    end

    return q0_out
end

"""
    checkSTM(dynamicsModel; relTol)

Return true if STM is accurate

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `relTol::Float64`: Relative tolerance (default = 2E-3)
"""
function checkSTM(dynamicsModel::BCR4BP41DynamicsModel, relTol::Float64 = 2E-3)
    stepSize::Float64 = sqrt(eps(Float64))
    numStates::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    propagator = MBD.Propagator()
    propagatorSTM = MBD.Propagator(equationType = MBD.STM)
    X::Vector{Float64} = [0.9, 0, 0, 0, -0.7, 0, 0]
    tau::Float64 = 0.1
    arc::MBD.BCR4BP41Arc = propagate(propagatorSTM, appendExtraInitialConditions(dynamicsModel, X, MBD.STM), [0, tau], dynamicsModel)
    STMAnalytical::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64}(getStateTransitionMatrix(dynamicsModel, getStateByIndex(arc, -1)))
    STMNumerical::StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64}(zeros(Float64, (numStates, numStates)))
    for index::Int16 in Int16(1):numStates
        perturbedFreeVariables::Vector{Float64} = copy(X)
        perturbedFreeVariables[index] -= stepSize
        arcMinus::MBD.BCR4BP41Arc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorMinus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcMinus, -1))
        perturbedFreeVariables[index] += 2*stepSize
        arcPlus::MBD.BCR4BP41Arc = propagate(propagator, perturbedFreeVariables, [0, tau], dynamicsModel)
        constraintVectorPlus::StaticArrays.SVector{Int64(numStates), Float64} = StaticArrays.SVector{Int64(numStates), Float64}(getStateByIndex(arcPlus, -1))
        STMNumerical[:,index] = (constraintVectorPlus-constraintVectorMinus)./(2*stepSize)
    end
    absDiff::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64}(STMNumerical-STMAnalytical)
    relDiff::StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64} = StaticArrays.MMatrix{Int64(numStates), Int64(numStates), Float64}(copy(absDiff))
    for r::Int16 in Int16(1):numStates
        for c::Int16 in Int16(1):numStates
            if abs(STMAnalytical[r,c]) < stepSize*1E3
                relDiff[r,c] = absDiff[r,c]
            elseif abs(STMNumerical[r,c]) > 1E-11
                relDiff[r,c] = absDiff[r,c]/STMNumerical[r,c]
            end
            if abs(relDiff[r,c]) > relTol
                throw(ErrorException("Jacobian error in entry ($r, $c): Expected = $(STMNumerical[r,c]); Actual = $(STMAnalytical[r,c]) (Relative error = $(relDiff[r,c]))"))
                return false
            end
        end
    end

    return true
end

"""
    evaluateEquations(dynamicsModel, equationType, t, q)

Return time derivative of state vector

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `equationType::EquationType`: EOM type
- `t::Float64`: Time [ndim]
- `q::Vector{Float64}`: State vector [ndim]
"""
function evaluateEquations(dynamicsModel::BCR4BP41DynamicsModel, equationType::MBD.EquationType, t::Float64, q::Vector{Float64})
    qdot::Vector{Float64} = Vector{Float64}(undef, getStateSize(dynamicsModel, equationType))
    EOMs::MBD.BCR4BP41EquationsOfMotion = getEquationsOfMotion(dynamicsModel, equationType)
    computeDerivatives!(qdot, q, (EOMs,), t)

    return qdot
end

"""
    getEpochDependencies(dynamicsModel, q)

Return derivative of state with respect to epoch

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getEpochDependencies(dynamicsModel::BCR4BP41DynamicsModel, q_full::Vector{Float64})
    n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
    (Int16(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)

    isEpochIndependent(dynamicsModel) ? (return zeros(Float64, n_simple)) : (return q_full[n_simple^2+1:n_simple^2+n_simple])
end

"""
    getEpochTime(dynamicsModel, frame, initialEpochGuess, theta20)

Return next corresponding epoch time

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `frame::String`: Fixed ecliptic frame
- `initialEpochGuess::String`: Initial epoch guess
- `theta40::Float64`: P2 angle [ndim]
"""
function getEpochTime(dynamicsModel::BCR4BP41DynamicsModel, frame::String, initialEpochGuess::String, theta20::Float64)
    tstar41::Float64 = get41CharTime(dynamicsModel)
    theta2dot::Float64 = evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.9, 0, 0, 0, -0.3, 0, 0])[7]
    epochTimeGuess::Float64 = SPICE.str2et(initialEpochGuess)
    theta2Diff::Float64 = pi
    while abs(theta2Diff) > 1E-8
        Q2::Vector{Float64} = getEphemerides(initialEpochGuess, [0.0], dynamicsModel.systemData.primaryNames[2], dynamicsModel.systemData.primaryNames[4], frame)[1][1]
        Q4::Vector{Float64} = getEphemerides(initialEpochGuess, [0.0], dynamicsModel.systemData.primaryNames[3], dynamicsModel.systemData.primaryNames[4], frame)[1][1]
        B1::MBD.BodyData = dynamicsModel.systemData.primaryData[4]
        P2SPICEElements::StaticArrays.MVector{20, Float64} = StaticArrays.MVector{20, Float64}(SPICE.oscltx(Q2, epochTimeGuess, B1.gravParam))
        P2SPICEElements[3] = 0.0
        R2::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(SPICE.conics(P2SPICEElements[1:8], epochTimeGuess)[1:3])
        P4SPICEElements::StaticArrays.MVector{20, Float64} = StaticArrays.MVector{20, Float64}(SPICE.oscltx(Q4, epochTimeGuess, B1.gravParam))
        P4SPICEElements[3] = 0.0
        R4::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(SPICE.conics(P4SPICEElements[1:8], epochTimeGuess)[1:3])
        r2::Float64 = LinearAlgebra.norm(R2)
        r4::Float64 = LinearAlgebra.norm(R4)
        theta2Guess::Float64 = pi-acos(LinearAlgebra.dot(R2, R4)/r2/r4)
        theta2Diff = theta20-theta2Guess
        tDiff::Float64 = (theta2Diff/theta2dot)*tstar41
        epochTimeGuess += tDiff
        initialEpochGuess = SPICE.et2utc(epochTimeGuess, :C, 11)
    end
    
    return epochTimeGuess
end

"""
    getEquationsOfMotion(dynamicsModel, equationType)

Return EOMs

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `equationType::EquationType`: EOM type
"""
function getEquationsOfMotion(dynamicsModel::BCR4BP41DynamicsModel, equationType::MBD.EquationType)
    return MBD.BCR4BP41EquationsOfMotion(equationType, dynamicsModel)
end

"""
    getExcursion(dynamicsModel, primary, q)

Return distance from primary

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `primary::Int64`: Primary identifier
- `q::Vector{Float64}`: State vector [ndim]
"""
function getExcursion(dynamicsModel::BCR4BP41DynamicsModel, primary::Int64, q::Vector{Float64})
    lstar41::Float64 = get41CharLength(dynamicsModel)
    primaryPos::Vector{Float64} = getPrimaryState(dynamicsModel, primary)[1:3]

    return LinearAlgebra.norm(q[1:3]-primaryPos)*lstar41
end

"""
    getHamiltonian(dynamicsModel, q)

Return BCR4BP P4-B1 Hamiltonian

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q::Vector{Float64}`: State vector [ndim]
"""
function getHamiltonian(dynamicsModel::BCR4BP41DynamicsModel, q::Vector{Float64})
    mu12::Float64 = get12MassRatio(dynamicsModel)
    mu41::Float64 = get41MassRatio(dynamicsModel)
    omm12::Float64 = 1-mu12
    omm41::Float64 = 1-mu41
    a4::Float64 = get4Distance(dynamicsModel)
    v_2::Float64 = q[4]^2+q[5]^2+q[6]^2
    x1::Float64 = omm41-mu12*cos(q[7])/a4
    y1::Float64 = -mu12*sin(q[7])/a4
    x2::Float64 = omm41+omm12*cos(q[7])/a4
    y2::Float64 = omm12*sin(q[7])/a4
    r_13::Float64 = sqrt((q[1]-x1)^2+(q[2]-y1)^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-x2)^2+(q[2]-y2)^2+q[3]^2)
    r_43::Float64 = sqrt((q[1]+mu41)^2+q[2]^2+q[3]^2)
    U::Float64 = mu14*omm12/r_13+mu41*mu12/r_23+omm41/r_43+(1/2)*(q[1]^2+q[2]^2)

    return 2*U-v_2
end

"""
    getParameterDependencies(dynamicsModel, q_full)

Return derivative of state with respect to parameters

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q_full::Vector{Float64}`: Full state vector [ndim]
"""
function getParameterDependencies(dynamicsModel::BCR4BP41DynamicsModel, q_full::Vector{Float64})
    n_full::Int16 = getStateSize(dynamicsModel, MBD.FULL)
    (Int16(length(q_full)) < n_full) && throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    i0::Int16 = n_simple*(n_simple+1)
    n::Int16 = n_full-i0
    (n == Int16(0)) && (return zeros(Float64, (n_simple,0)))
    n_params::Int16 = n/n_simple
    dqdp::Matrix{Float64} = zeros(Float64, (n_simple,n_params))
    for r::Int16 in Int16(1):n_simple
        for c::Int16 in Int16(1):n_params
            dqdp[r,c] = q_full[i0+n_simple*(c-1)+r]
        end
    end

    return dqdp
end

"""
    getPrimaryState(dynamicsModel, primary, theta2)

Return state of primary in rotating frame

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `primary::Int64`: Primary identifier
- `theta2::Float64`: P2 angle [ndim]
"""
function getPrimaryState(dynamicsModel::BCR4BP41DynamicsModel, primary::Int64, theta2::Float64)
    (1 <= primary <= 2) || (primary == 4) || throw(ArgumentError("Invalid primary $primary"))
    mu12::Float64 = get12MassRatio(dynamicsModel)
    mu41::Float64 = get41MassRatio(dynamicsModel)
    omm12::Float64 = 1-mu12
    a4::Float64 = get4Distance(dynamicsModel)
    theta2dot::Float64 = evaluateEquations(dynamicsModel, MBD.SIMPLE, 0.0, [0.9, 0, 0, 0, -0.3, 0, 0])[7]
    q::Vector{Float64} = push!(zeros(Float64, 6), theta2)
    if primary == 1
        q[1] = 1-mu41-mu12*cos(q[7])/a4
        q[2] = -mu12*sin(q[7])/a4
        q[4] = mu12*theta2dot*sin(q[7])/a4
        q[5] = -mu12*theta2dot*cos(q[7])/a4
    elseif primary == 2
        q[1] = 1-mu41+omm12*cos(q[7])/a4
        q[2] = omm12*sin(q[7])/a4
        q[4] = -omm12*theta2dot*sin(q[7])/a4
        q[5] = omm12*theta2dot*cos(q[7])/a4
    elseif primary == 4
        q[1] = -mu41
    end

    return q
end

"""
    getPsuedopotentialJacobian(dynamicsModel, q)

Return second derivative of pseudopotential function at given location and P2 angle

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q::Vector{Float64}`: State vector [ndim]
"""
function getPseudopotentialJacobian(dynamicsModel::BCR4BP41DynamicsModel, q::Vector{Float64})
    mu12::Float64 = get12MassRatio(params[1])
    mu41::Float64 = get41MassRatio(params[1])
    omm12::Float64 = 1-mu12
    omm41::Float64 = 1-mu41
    a4::Float64 = get4Distance(params[1])
    x1::Float64 = omm41-mu12*cos(q[7])/a4
    y1::Float64 = -mu12*sin(q[7])/a4
    x2::Float64 = omm41+omm12*cos(q[7])/a4
    y2::Float64 = omm12*sin(q[7])/a4
    r_13::Float64 = sqrt((q[1]-x1)^2+(q[2]-y1)^2+q[3]^2)
    r_23::Float64 = sqrt((q[1]-x2)^2+(q[2]-y2)^2+q[3]^2)
    r_43::Float64 = sqrt((q[1]+mu41)^2+q[2]^2+q[3]^2)
    r_13_3::Float64 = r_13^3
    r_23_3::Float64 = r_23^3
    r_43_3::Float64 = r_43^3
    ddUdr::Vector{Float64} = zeros(Float64, 9)
    ddUdr[1] = 1-mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*(q[1]-x1)^2/r_13_5+3*mu41*mu12*(q[1]-x2)^2/r_23_5+3*omm41*(q[1]+mu41)^2/r_43_5
    ddUdr[2] = 1-mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*(q[2]-y1)^2/r_13_5+3*mu41*mu12*(q[2]-y2)^2/r_23_5+3*omm41*q[2]^2/r_43_5
    ddUdr[3] = -mu41*omm12/r_13_3-mu14*mu12/r_23_3-omm41/r_43_3+3*mu41*omm12*q[3]^2/r_13_5+3*mu41*mu12*q[3]^2/r_23_5+3*omm41*q[3]^2/r_43_5
    ddUdr[4] = 3*mu41*omm12*(q[1]-x1)*(q[2]-y1)/r_13_5+3*mu41*mu12*(q[1]-x2)*(q[2]-y2)/r_23_5+3*omm41*(q[1]+mu41)*q[2]/r_43_5
    ddUdr[5] = 3*mu41*omm12*(q[1]-x1)*q[3]/r_13_5+3*mu41*mu12*(q[1]-x2)*q[3]/r_23_5+3*omm41*(q[1]+mu41)*q[3]/r_43_5
    ddUdr[6] = 3*mu41*omm12*(q[2]-y1)*q[3]/r_13_5+3*mu41*mu12*(q[2]-y2)*q[3]/r_23_5+3*omm41*q[2]*q[3]/r_43_5
    ddUdr[7] = mu41*mu12*omm12*sin(q[7])/(a4*r_13_3)-mu41*mu12*omm12*sin(q[7])/(a4*r_23_3)-3*mu41*mu12*omm12*(q[1]-x1)*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu14*mu12*omm12*(q[1]-x2)*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)
    ddUdr[8] = -mu41*mu12*omm12*cos(q[7])/(a4*r_13_3)+mu41*mu12*omm12*cos(q[7])/(a4*r_23_3)-3*mu41*mu12*omm12*(q[2]-y1)*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu41*mu12*omm12*(q[2]-y2)*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)
    ddUdr[9] = -3*mu41*mu12*omm12*q[3]*((q[1]-x1)*sin(q[7])-(q[2]-y1)*cos(q[7]))/(a4*r_13_5)+3*mu41*mu12*omm12*q[3]*((q[1]-x2)*sin(q[7])-(q[2]-y2)*cos(q[7]))/(a4*r_23_5)

    return ddUdr
end

"""
    getStateSize(dynamicsModel, equationType)

Return number of state variables

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `equationType::EquationType`: EOM type
"""
function getStateSize(dynamicsModel::BCR4BP41DynamicsModel, equationType::MBD.EquationType)
    type = Dict(MBD.SIMPLE => Int16(7), MBD.STM => Int16(56), MBD.FULL => Int16(56), MBD.ARCLENGTH => Int16(57), MBD.MOMENTUM => Int16(57))

    return type[equationType]
end

"""
    getStateTransitionMatrix(dynamicsModel, q0)

Return STM

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `q0::Vector{Float64}`: Initial state vector with STM in column-major order [ndim]
"""
function getStateTransitionMatrix(dynamicsModel::BCR4BP41DynamicsModel, q0::Vector{Float64})
    n_STM::Int16 = getStateSize(dynamicsModel, MBD.STM)
    (Int16(length(q0)) < n_STM) && throw(ArgumentError("State vector length is $(length(q0)), but should be at least $n_STM"))
    n_simple::Int16 = getStateSize(dynamicsModel, MBD.SIMPLE)
    STM::Matrix{Float64} = zeros(Float64, (n_simple,n_simple))
    for r::Int16 in Int16(1):n_simple
        for c::Int16 in Int16(1):n_simple
            STM[r,c] = q0[n_simple+n_simple*(c-1)+r]
        end
    end

    return STM
end

"""
    gettheta2(dynamicsModel, frame, initialEpoch)

Return P2 angle [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `frame::String`: Fixed ecliptic frame
- `initialEpoch::String`: Initial epoch
"""
function gettheta2(dynamicsModel::BCR4BP41DynamicsModel, frame::String, initialEpoch::String)
    epochTime::Float64 = SPICE.str2et(initialEpochGuess)
    Q2::Vector{Float64} = getEphemerides(initialEpoch, [0.0], dynamicsModel.systemData.primaryNames[2], dynamicsModel.systemData.primaryNames[4], frame)[1][1]
    Q4::Vector{Float64} = getEphemerides(initialEpoch, [0.0], dynamicsModel.systemData.primaryNames[3], dynamicsModel.systemData.primaryNames[4], frame)[1][1]
    B1::MBD.BodyData = dynamicsModel.systemData.primaryData[4]
    P2SPICEElements::StaticArrays.MVector{20, Float64} = StaticArrays.MVector{20, Float64}(SPICE.oscltx(Q2, epochTime, B1.gravParam))
    P2SPICEElements[3] = 0.0
    R2::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(SPICE.conics(P2SPICEElements[1:8], epochTime)[1:3])
    P4SPICEElements::StaticArrays.MVector{20, Float64} = StaticArrays.MVector{20, Float64}(SPICE.oscltx(Q4, epochTime, B1.gravParam))
    P4SPICEElements[3] = 0.0
    R4::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(SPICE.conics(P4SPICEElements[1:8], epochTime)[1:3])
    r2::Float64 = LinearAlgebra.norm(R2)
    r4::Float64 = LinearAlgebra.norm(R4)

    return pi-acos(LinearAlgebra.dot(R2, R4)/r2/r4)
end

"""
    get12MassRatio(dynamicsModel)

Return BCR4BP P1-P2 system mass ratio

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get12MassRatio(dynamicsModel::BCR4BP41DynamicsModel)
    return get12MassRatio(dynamicsModel.systemData)
end

"""
    get2BApproximation(dynamicsModel, frame, primary, radius, initialEpochGuess, theta2)

Return states of 2BP approximation about primary [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `frame::String`: Fixed ecliptic frame
- `primary::Int64`: Primary identifier
- `radius::Float64`: Circular radius [ndim]
- `initialEpochGuess::String`: Initial epoch guess
- `theta2::Float64`: Desired P2 angle [ndim]
"""
function get2BApproximation(dynamicsModel::BCR4BP41DynamicsModel, frame::String, primary::Int64, radius::Float64, initialEpochGuess::String, theta2::Float64)
    lstar41::Float64 = get41CharLength(dynamicsModel)
    tstar41::Float64 = get41CharTime(dynamicsModel)
    radiusDim::Float64 = radius*lstar41
    vDim::Float64 = sqrt(bodyData.gravParam/radiusDim)
    v::Float64 = vDim*tstar41/lstar41
    q_primaryInertial::Vector{Float64} = [-radius, 0, 0, 0, v, 0]
    t::Float64 = getEpochTime(dynamicsModel, frame, initialEpochGuess, theta2)

    return primaryEclipticToRotating41(dynamicsModel, frame, primary, [q_primaryInertial], [t])[1][1]
end

"""
    get4Distance(dynamicsModel)

Return BCR4BP P4 distance from B1 [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get4Distance(dynamicsModel::BCR4BP41DynamicsModel)
    return get4Distance(dynamicsModel.systemData)
end

"""
    get4Mass(dynamicsModel)

Return BCR4BP P4 mass [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get4Mass(dynamicsModel::BCR4BP41DynamicsModel)
    return get4Mass(dynamicsModel.systemData)
end

"""
    get41CharLength(dynamicsModel)

Return BCR4BP P4-B1 characteristic length

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get41CharLength(dynamicsModel::BCR4BP41DynamicsModel)
    return get41CharLength(dynamicsModel.systemData)
end

"""
    get41CharTime(dynamicsModel)

Return BCR4BP P4-B1 characteristic time

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get41CharTime(dynamicsModel::BCR4BP41DynamicsModel)
    return get41CharTime(dynamicsModel.systemData)
end

"""
    get41MassRatio(dynamicsModel)

Return BCR4BP P4-B1 system mass ratio

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function get41MassRatio(dynamicsModel::BCR4BP41DynamicsModel)
    return get41MassRatio(dynamicsModel.systemData)
end

"""
    isEpochIndependent(dynamicsModel)

Return true if dynamics model is epoch independent

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
function isEpochIndependent(dynamicsModel::BCR4BP41DynamicsModel)
    return false
end

"""
    primaryEclipticToRotating41(dynamicsModel, frame, primary, states, times)

Return BCR4BP P4-B1 rotating frame states and times [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `frame::String`: Fixed ecliptic frame
- `primary::Int64`: Primary identifier
- `states::Vector{Vector{Float64}}`: Primary-centered fixed states [ndim]
- `times::Vector{Float64}`: Epoch times [s]
"""
function primaryEclipticToRotating41(dynamicsModel::BCR4BP41DynamicsModel, frame::String, primary::Int64, states::Vector{Vector{Float64}}, times::Vector{Float64})
    dynamicsModel12 = MBD.BCR4BP12DynamicsModel(dynamicsModel.systemData)
    (states12::Vector{Vector{Float64}}, times12::Vector{Float64}) = primaryEclipticToRotating12(dynamicsModel12, frame, primary, states, times)

    return rotating12ToRotating41(dynamicsModel12, states12, times12)
end

"""
    rotating41ToPrimaryEcliptic(dynamicsModel, frame, primary, initialEpochGuess, states, times)

Return primary-centered fixed frame states and times [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `frame::String`: Fixed ecliptic frame
- `primary::Int64`: Primary identifier
- `initialEpochGuess::String`: Initial epoch guess
- `states::Vector{Vector{Float64}}`: Rotating states [ndim]
- `times::Vector{Float64}`: Rotating times [ndim]
"""
function rotating41ToPrimaryEcliptic(dynamicsModel::BCR4BP41DynamicsModel, frame::String, primary::Int64, initialEpochGuess::String, states::Vector{Vector{Float64}}, times::Vector{Float64})
    (states12::Vector{Vector{Float64}}, times12::Vector{Float64}) = rotating41ToRotating12(dynamicsModel, states, times)
    dynamicsModel12 = MBD.BCR4BP12DynamicsModel(dynamicsModel.systemData)
    
    return rotating12ToPrimaryEcliptic(dynamicsModel12, frame, primary, initialEpochGuess, states12, times12)
end

"""
    rotating41ToRotating12(dynamicsModel, states41, times41)

Return BCR4BP P1-P2 rotating frame states and times [ndim]

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `states41::Vector{Vector{Float64}}`: BCR4BP P4-B1 rotating frame states [ndim]
- `times41::Vector{Float64}`: BCR4BP P4-B1 rotating frame times [ndim]
"""
function rotating41ToRotating12(dynamicsModel::BCR4BP41DynamicsModel, states41::Vector{Vector{Float64}}, times41::Vector{Float64})
    numTimes::Int16 = Int16(length(times41))
    m4::Float64 = get4Mass(dynamicsModel.systemData)
    a4::Float64 = get4Distance(dynamicsModel.systemData)
    theta2dot::Float64 = sqrt((a4^3)/(m4+1))-1
    states12::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
    for t::Int16 = Int16(1):numTimes
        state::StaticArrays.SVector{7, Float64} = StaticArrays.SVector{7, Float64}(states41[t])
        theta2::Float64 = state[7]
        C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([cos(theta2) sin(theta2) 0; -sin(theta2) cos(theta2) 0; 0 0 1])
        Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}(theta2dot.*[-sin(theta2) cos(theta2) 0; -cos(theta2) -sin(theta2) 0; 0 0 0])
        states12[t] = [a4.*C zeros(Float64, 3, 4); sqrt((m4+1)/a4).*Cdot sqrt((m4+1)/a4).*C zeros(Float64, 3, 1); zeros(Float64, 1, 6) -1]*(state+append!([-1+get41MassRatio(dynamicsModel.systemData)], zeros(Float64, 6)))+append!(zeros(Float64, 6), [pi])
    end
    times12::Vector{Float64} = (get41CharTime(dynamicsModel.systemData)/get12CharTime(dynamicsModel.systemData)).*times41

    return (states12, times12)
end
