"""
DynamicsModel methods

Author: Jonathan LeFevre Richmond
C: 5/4/26
U: 6/26/26

QUEUE:
    checkSTM() needs Propagator, propagate(), getStateTransitionMatrix(), Arc, getStateByIndex()
    evaluateEquations() needs EquationsOfMotion, getEquationsOfMotion(), computeDerivatives!()
    getEpochDependencies() needs isEpochIndependent()
    getEquationsOfMotion() needs EquationsOfMotion
    getLinearVariation() needs getPseudopotentialJacobian

TO DO:
"""


"""
    adjustInitialConditions(dynamicsModel::AbstractDynamicsModel, q0::AbstractVector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType) -> Vector{Float64}

Return initial conditions for output equations of motion type

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q0::AbstractVector{Float64}`: Initial conditions
- `inputEquationType::EquationType`: Equations of motion type for `q0`
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector{Float64}`: Initial conditions

Errors
- Throws `MethodError` if not implemented for the dynamics model

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, the only non-zero values are in the simple state
    and STM

Example
```
q0_full::Vector{Float64} = adjustInitialConditions(dynamicsModel, q0_STM, STM, FULL)
```
"""
function adjustInitialConditions(dynamicsModel::AbstractDynamicsModel, q0::AbstractVector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered generic adjustInitialConditions" dynamicsModel inputEquationType outputEquationType

    Logging.@error "adjustInitialConditions is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(adjustInitialConditions, (dynamicsModel, q0, inputEquationType, outputEquationType)))
end

"""
    appendExtraInitialConditions(dynamicsModel::AbstractDynamicsModel, q0_simple::AbstractVector{Float64}, outputEquationType::EquationType) -> Vector{Float64}

Return initial conditions for output equations of motion type

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q0_simple::AbstractVector{Float64}`: Simple initial conditions
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector{Float64}`: Initial conditions

Errors
- No error checking

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when forwarding to
    `adjustInitialConditions()`

Notes
- This is a common special case of `adjustInitialConditions` and forwards to
    that method
- For `CR3BPDynamicsModel`, the only non-zero values are in the simple state
    and STM

Example
```
q0_full::Vector{Float64} = appendExtraInitialConditions(dynamicsModel, q0_simple, FULL)
```
"""
function appendExtraInitialConditions(dynamicsModel::AbstractDynamicsModel, q0_simple::AbstractVector{Float64}, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered generic appendExtraInitialConditions" dynamicsModel outputEquationType

    Logging.@debug "Forwarding to adjustInitialConditions" dynamicsModel outputEquationType
    return adjustInitialConditions(dynamicsModel, q0_simple, SIMPLE, outputEquationType)
end

"""
    extractStateTransitionMatrix(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64}) -> Matrix{Float64}

Return STM

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q::AbstractVector{Float64}`: State vector

Returns
- `Matrix{Float64}`: State transition matrix

Errors
- Throws `ArgumentError` if state vector is non-finite or too short
- Throws `DomainError` if STM Frobenius norm is not finite

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when extracting STM, or when
    STM is returned

Example
```
Φ::matrix{Float64} = extractStateTransitionMatrix(dynamicsModel, q)
```
"""
function extractStateTransitionMatrix(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64})::Matrix{Float64}
    Logging.@debug "Entered extractStateTransitionMatrix" dynamicsModel q

    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    n_in::Int64 = length(q)
    n_STM::Int64 = getStateSize(dynamicsModel, STM)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)

    # Validate that state vector is long enough to contain full STM block
    if n_in < n_STM
        Logging.@error "State vector is too short to contain STM block" required=n_STM actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_STM to extract STM"))
    end

    # Extract and reshape flattened STM block into matrix
    Φ::Matrix{Float64} = reshape(q[n_simple+1:n_STM], n_simple, n_simple)
        
    Logging.@debug "Returning STM" frobeniusNorm=frobenius

    return Φ
end

"""
    getCharLengths(dynamicsModel::AbstractDynamicsModel)

Return dynamics model characteristic length scales

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object

Returns
- Variable

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns secondary body's orbital radius [km]

Example
```
lstars = getCharLengths(dynamicsModel)
```
"""
function getCharLengths(dynamicsModel::AbstractDynamicsModel)
    Logging.@debug "Entered generic getCharLengths" dynamicsModel

    Logging.@error "getCharLengths is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getCharLengths, (dynamicsModel,)))
end

"""
    getCharMasses(dynamicsModel::AbstractDynamicsModel)

Return dynamics model characteristic mass scales

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object

Returns
- Variable

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns sum of gravitational parameters divided by
    the gravitational constant [kg]

Example
```
mstars = getCharMasses(dynamicsModel)
```
"""
function getCharMasses(dynamicsModel::AbstractDynamicsModel)
    Logging.@debug "Entered generic getCharMasses" dynamicsModel

    Logging.@error "getCharMasses is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getCharMasses, (dynamicsModel,)))
end

"""
    getCharTimes(dynamicsModel::AbstractDynamicsModel)

Return dynamics model characteristic time scales

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object

Returns
- Variable

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns orbital period scale derived from Kepler's
    3rd law [s]

Example
```
tstars = getCharTimes(dynamicsModel)
```
"""
function getCharTimes(dynamicsModel::AbstractDynamicsModel)
    Logging.@debug "Entered generic getCharTimes" dynamicsModel

    Logging.@error "getCharTimes is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getCharTimes, (dynamicsModel,)))
end

"""
    getEquilibriumPoint(dynamicsModel::AbstractDynamicsModel, point::Int64) -> Vector{Float64}

Return equilibrium point state

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `point::Int64`: Equilibrium point index

Returns
- `Vector{Float64}`: Equilibrium point state

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns Lagrange point state [ndim]

Example
```
q_L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
```
"""
function getEquilibriumPoint(dynamicsModel::AbstractDynamicsModel, point::Int64)::Vector{Float64}
    Logging.@debug "Entered generic getEquilibrumPoint" dynamicsModel

    Logging.@error "getEquilibriumPoint is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getEquilibriumPoint, (dynamicsModel, point)))
end

"""
    getExcursion(dynamicsModel::AbstractDynamicsModel, primary::Int64, q::AbstractVector{Float64}) -> Float64

Return excursion from primary

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `primary::Int64`: Primary index
- `q::AbstractVector{Float64}`: State vector

Returns
- `Float64`: Excursion distance

Errors
- Throws `ArgumentError` if primary index is invlaid or if state vector is
    non-finite or too short

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when excursion is returned

Example
```
d::Float64 = getExcursion(dynamicsModel, 1, q)
```
"""
function getExcursion(dynamicsModel::AbstractDynamicsModel, primary::Int64, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered getExcursion" dynamicsModel primary q

    # Validate primary index
    if !(1 <= primary <= getNumPrimaries(dynamicsModel))
        Logging.@error "Invalid primary index" primary=primary nPrimaries=getNumPrimaries(dynamicsModel)
        throw(ArgumentError("Primary must be between 1 and $(getNumPrimaries(dynamicsModel)), got $primary"))
    end
    
    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end
  
    # Validate that q contains at least position components
    if length(q) < 3
        Logging.@error "State vector is too short to contain position components" rquired=3 actual=length(q)
        throw(ArgumentError("State vector length $(length(q)) is insufficient; need at least 3"))
    end

    q_P::Vector{Float64} = getPrimaryState(dynamicsModel, primary)

    # Compute Euclidean distance directly
    dx::Float64 = q[1]-q_P[1]
    dy::Float64 = q[2]-q_P[2]
    dz::Float64 = q[3]-q_P[3]
    d::Float64 = sqrt(dx^2+dy^2+dz^2)

    Logging.@debug "Returning excursion" primary=primary excursion=d

    return d
end

"""
    getHamiltonian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64}) -> Float64

Return Hamiltonian value

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q::AbstractVector{Float64}`: State vector

Returns
- `Float64`: Hamiltonian value

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns Jacobi constant [ndim]

Example
```
H::Float64 = getHamiltonian(dynamicsModel, q)
```
"""
function getHamiltonian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered generic getHamiltonian" dynamicsModel

    Logging.@error "getHamiltonian is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getHamiltonian, (dynamicsModel, q)))
end

"""
    getMassRatios(dynamicsModel::AbstractDynamicsModel)

Return dynamics model mass ratios

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object

Returns
- Variable

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns ratio between secondary and total mass

Example
```
μs = getMassRatios(dynamicsModel)
```
"""
function getMassRatios(dynamicsModel::AbstractDynamicsModel)
    Logging.@debug "Entered generic getMassRatios" dynamicsModel

    Logging.@error "getMassRatios is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getMassRatios, (dynamicsModel,)))
end

"""
    getNumPrimaries(dynamicsModel::AbstractDynamicswModel) -> Int64

Return number of primary bodies contained in an `AbstractDynamicsModel` object

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object

Returns
- `Int64`: Number of entries in `dynamicsModel.primaryData`

Errors
- None

Logging
- Emits `@warn` logs if `dynamicsModel.primaryData` is empty
- Emits `@debug` logs when function is entered or when
    `dynamicsModel.primaryData` contains entries

Example
```
nPrimaries::Int64 = getNumPrimaries(dynamicsModel)
```
"""
function getNumPrimaries(dynamicsModel::AbstractDynamicsModel)::Int64
    Logging.@debug "Entered getNumPrimaries" dynamicsModel

    n::Int64 = length(dynamicsModel.primaryData)
    if n == 0
        Logging.@warn "AbstractDynamicsModel contains no primaries (primaryData is empty)"
    else
        Logging.@debug "Returning primary count" count=n
    end
    
    return n
end

"""
    getParameterDependencies(dynamicsModel::AbstractDynamicsModel, q_full::Vector{Float64}) -> Matrix{Float64}

Return derivative of state with respect to parameters

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q_full::Vector{Float64}`: Full state vector

Returns
- `Matrix{Float64}`: Derivative of state with respect to parameters

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns empty matrix

Example
```
dqdp::Matrix{Float64} = getParameterDependencies(dynamicsModel, q_full)
```
"""
function getParameterDependencies(dynamicsModel::AbstractDynamicsModel, q_full::Vector{Float64})::Matrix{Float64}
    Logging.@debug "Entered generic getParameterDependencies" dynamicsModel

    Logging.@error "getParameterDependencies is not implemented for this dynamicsModel type" type=typeof(dynamicsModel)
    throw(MethodError(getParameterDependencies, (dynamicsModel, q_full)))
end

"""
    getPrimaryState(dynamicsModel::AbstractDynamicsModel, primary::Int64) -> Vector{Float64}

Return primary state

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `primary::Int64`: Primary index

Returns
- `Vector{Float64}`: Primary state

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns primary state [ndim]

Example
```
q_1::Vector{Float64} = getPrimaryState(dynamicsModel, 1)
```
"""
function getPrimaryState(dynamicsModel::AbstractDynamicsModel, primary::Int64)::Vector{Float64}
    Logging.@debug "Entered generic getPrimaryState" dynamicsModel

    Logging.@error "getPrimaryState is not implemented for this dynamicsModel type" type=typeof(dynamicsModel)
    throw(MethodError(getPrimaryState, (dynamicsModel, primary)))
end

"""
    getPseudopotential(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64}) -> Float64

Return pseudo-potential

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q::AbstractVector{Float64}`: State vector

Returns
- `Float64`: Pseudo-potential

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns pseudo-potential [ndim]

Example
```
U::Float64 = getPseudopotential(dynamicsModel, q)
```
"""
function getPseudopotential(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered generic getPseudopotential" dynamicsModel

    Logging.@error "getPseudopotential is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getPseudopotential, (dynamicsModel, q)))
end

"""
    getPseudopotentialHessian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64}) -> Vector{Float64}

Return pseudo-potential Hessian

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q::AbstractVector{Float64}`: State vector

Returns
- `Vector{Float64}`: Pseudo-potential Hessian elements

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns pseudo-potential Hessian elements [ndim]

Example
```
d2Udr2::Float64 = getPseudopotentialHessian(dynamicsModel, q)
```
"""
function getPseudopotentialHessian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64})::Vector{Float64}
    Logging.@debug "Entered generic getPseudopotentialHessian" dynamicsModel

    Logging.@error "getPseudopotentialHessian is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getPseudopotentialHessian, (dynamicsModel, q)))
end

"""
    getPseudopotentialJacobian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64}) -> Vector{Float64}

Return pseudo-potential Jacobian

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q::AbstractVector{Float64}`: State vector

Returns
- `Vector{Float64}`: Pseudo-potential Jacobian

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type
- For `CR3BPDynamicsModel`, returns pseudo-potential jacobian [ndim]

Example
```
dUdr::Float64 = getPseudopotentialJacobian(dynamicsModel, q)
```
"""
function getPseudopotentialJacobian(dynamicsModel::AbstractDynamicsModel, q::AbstractVector{Float64})::Vector{Float64}
    Logging.@debug "Entered generic getPseudopotentialJacobian" dynamicsModel

    Logging.@error "getPseudopotentialJacobian is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getPseudopotentialJacobian, (dynamicsModel, q)))
end

"""
    getStateSize(dynamicsModel::AbstractDynamicsModel, equationType::EquationType) -> Int64

Return state vector size for equations of motion type

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `equationType::EquationType`: Equations of motion type

Returns
- `Int64`: State vector size

Errors
- Throws `MethodError` if not implemented for the dynamics model type

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Notes
- This is a generic method that dispatches on dynamics model type

Example
```
n_full::Int64 = getStateSize(dynamicsModel, FULL)
```
"""
function getStateSize(dynamicsModel::AbstractDynamicsModel, equationType::EquationType)::Int64
    Logging.@debug "Entered generic getStateSize" dynamicsModel equationType

    Logging.@error "getStateSize is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(getStateSize, (dynamicsModel, equationType)))
end


"""
    adjustInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0::AbstractVector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType) -> Vector{Float64}

Return initial conditions for CR3BP output equations of motion type

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q0::AbstractVector{Float64}`: Initial conditions [ndim]
- `inputEquationType::EquationType`: Equations of motion type for `q0`
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector{Float64}`: Initial conditions [ndim]

Errors
- Throws `ArgumentError` if length of `q0` does not match that expected for
    `inputEquationType`

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when truncating or extending
    input initial conditions, when constructing STM initial conditions, or when
    adjusted initial conditions are returned

Example
```
q0_full::Vector{Float64} = adjustInitialConditions(dynamicsModel, q0_STM, STM, FULL)
```
"""
function adjustInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0::AbstractVector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered appendExtraInitialConditions" dynamicsModel inputEquationType outputEquationType

    # Validate input state vector length against expected size for inputEquationType
    n_in::Int64 = getStateSize(dynamicsModel, inputEquationType)
    if length(q0) != n_in
        Logging.@error "Input state vector has incorrect length" expected=n_in actual=length(q0)
        throw(ArgumentError("State vector length is $(length(q0)), but should be $n_in"))
    end

    # Validate input state vector
    if !all(isfinite, q0)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q0)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    n_out::Int64 = getStateSize(dynamicsModel, outputEquationType)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    n_STM::Int64 = getStateSize(dynamicsModel, STM)

    # Direct truncation
    if outputEquationType == SIMPLE
        Logging.@debug "Truncating state vector directly" inputType=inputEquationType outputType=outputEquationType inputSize=n_in outputSize=n_out

        return q0[1:n_out]
    end
    q0_out::Vector{Float64} = zeros(Float64, n_out)

    if outputEquationType == FULL
        # FULL layout must be handled explicitly
        Logging.@debug "Constructing full state vector" inputSize=n_in outputSize=n_out

        # Copy simple states from leading portion of input
        q0_out[1:n_simple] .= q0[1:n_simple]

        if n_in >= n_STM
            # Input already contains full STM block - copy it directly
            q0_out[n_simple+1:n_STM] .= q0[n_simple+1:n_STM]
        else
            # Input has no STM block - initialize it as a flattened identity matrix
            for j in n_simple+1:n_simple+1:n_STM
                q0_out[j] = 1.0
            end
        end
    elseif n_in >= n_out
        # Output is no larger than input - strip to SIMPLE, then build up
        Logging.@debug "Truncating state vector indirectly" inputType=inputEquationType outputType=outputEquationType inputSize=n_in outputSize=n_out
        q0_out[1:n_simple] .= q0[1:n_simple]
    else
        # Output is larger than input but is not FULL - extend with zeros, then populate
        Logging.@debug "Extending state vector" inputSize=n_in outputSize=n_out

        # Copy simple states into leading portion of output vector
        q0_out[1:n_simple] .= q0[1:n_simple]

        # Initialize STM block as flattened identity matrix only when input doesn't already contain STM block but output does
        if (n_in < n_STM) && (n_out >= n_STM)
            Logging.@debug "Constructing STM initial conditions" inputSize=n_in outputSize=n_out
            for j in n_simple+1:n_simple+1:n_STM
                q0_out[j] = 1.0
            end
        end
    end

    Logging.@debug "Returning adjusted state vector" outputSize=n_out
        
    return q0_out
end

"""
    getCharLengths(dynamicsModel::CR3BPDynamicsModel) -> Float64
    
Return CR3BP characteristic length scale

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object

Returns
- `Float64`: Secondary body's orbital radius [km]

Errors
- Throws `DomainError` if characteristic length is not finite or non-positive

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when characteristic length is
    returned

Example
```
lstar::Float64 = getCharLengths(dynamicsModel)
```
"""
function getCharLengths(dynamicsModel::CR3BPDynamicsModel)::Float64
    Logging.@debug "Entered getCharLengths (CR3BP)" dynamicsModel

    # Extract orbital radius for secondary body
    lstar::Float64 = dynamicsModel.primaryData[2].a
    if !isfinite(lstar)
        Logging.@error "Characteristic length is not finite" charLength=lstar
        throw(DomainError(lstar, "Characteristic length must be finite"))
    end
    if lstar <= 0.0
        Logging.@error "Characteristic length is non-positive" charLength=lstar
        throw(DomainError(lstar, "Characteristic length must be positive"))
    end

    Logging.@debug "Returning characteristic length" charLength=lstar

    return lstar
end

"""
    getCharMasses(dynamicsModel::CR3BPDynamicsModel) -> Float64
    
Return CR3BP characteristic mass scale

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object

Returns
- `Float64`: Sum of gravitational parameters divided by the gravitational
    constant [kg]

Errors
- Throws `DomainError` if gravitational parameters are not finite or
    non-positive

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when characteristic mass is
    returned

Example
```
mstar::Float64 = getCharMasses(dynamicsModel)
```
"""
function getCharMasses(dynamicsModel::CR3BPDynamicsModel)::Float64
    Logging.@debug "Entered getCharMasses (CR3BP)" dynamicsModel

    # Extract gravitational parameters for both primaries
    μ_1::Float64 = dynamicsModel.primaryData[1].μ
    μ_2::Float64 = dynamicsModel.primaryData[2].μ
    for (label, μ) in (("μ_1", μ_1), ("μ_2", μ_2))
        if !isfinite(μ)
            Logging.@error "Gravitational parameter is not finite" param=label value=μ
            throw(DomainError(μ, "$label must be finite"))
        end
        if μ <= 0.0
            Logging.@error "Gravitational parameter is non-positive" param=label value=μ
            throw(DomainError(μ, "$label must be positive"))
        end
    end

    # Calculate characteristic mass
    mstar::Float64 = (μ_1+μ_2)/GRAVITY

    Logging.@debug "Returning characteristic mass" charMass=mstar

    return mstar
end

"""
    getCharTimes(dynamicsModel::CR3BPDynamicsModel) -> Float64
    
Return CR3BP characteristic time scale

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object

Returns
- `Float64`: Orbital period scale derived from Kepler's 3rd law [s]

Errors
- Throws `DomainError` if gravitational parameters are not finite or
    non-positive

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when characteristic time is
    returned

Example
```
tstar::Float64 = getCharTimes(dynamicsModel)
```
"""
function getCharTimes(dynamicsModel::CR3BPDynamicsModel)::Float64
    Logging.@debug "Entered getCharTimes (CR3BP)" dynamicsModel

    # Extract characteristic length and gravitational parameters
    lstar::Float64 = getCharLengths(dynamicsModel)
    μ_1::Float64 = dynamicsModel.primaryData[1].μ
    μ_2::Float64 = dynamicsModel.primaryData[2].μ
    for (label, μ) in (("μ_1", μ_1), ("μ_2", μ_2))
        if !isfinite(μ)
            Logging.@error "Gravitational parameter is not finite" param=label value=μ
            throw(DomainError(μ, "$label must be finite"))
        end
        if μ <= 0.0
            Logging.@error "Gravitational parameter is non-positive" param=label value=μ
            throw(DomainError(μ, "$label must be positive"))
        end
    end

    # Calculate characteristic time
    tstar::Float64 = sqrt(lstar^3/(μ_1+μ_2))

    Logging.@debug "Returning characteristic time" charTime=tstar

    return tstar
end

"""
    getEquilibriumPoint(dynamicsModel::CR3BPDynamicsModel, point::Int64) -> Vector{Float64}

Return CR3BP equilibrium point state

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `point::Int64`: Equilibrium point index

Returns
- `Vector{Float64}`: Equilibrium point state [ndim]

Errors
- Throws `ArgumentError` if `point` is not between 1 and 5
- Throws `ErrorException` if Newton-Raphson algorithm does not converge within
    20 iterations
- Throws `DomainError` if equilibrium point state has non-finite elements

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing equilibrium
    point state, if using Newton-Raphson iterative process, or when equilibrium
    point state is returned

Example
```
q_L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
```
"""
function getEquilibriumPoint(dynamicsModel::CR3BPDynamicsModel, point::Int64)::Vector{Float64}
    Logging.@debug "Entered getEquilibriumPoint (CR3BP)" dynamicsModel point

    # Validate equilibrium point index
    if !(1 <= point <= 5)
        Logging.@error "Invalid equilibrium point index" point
        throw(ArgumentError("Equilibrium point must be between 1 and 5, got $point"))
    end

    μ::Float64 = getMassRatios(dynamicsModel)

    Logging.@debug "Computing equilibrium point state" point μ

    tol::Float64 = 1E-14
    maxCount::Int64 = 20
    q::Vector{Float64} = zeros(Float64, 6)

    # Newton-Raphson state
    γ::Float64 = 0.0
    γ_prev::Float64 = Inf
    count::Int64 = 0

    if point == 1
        Logging.@debug "Iterating on L1 position" μ

        # L1 point between two primaries
        # Initial guess from Hill's sphere approximation
        γ = cbrt(μ/(3(1-μ)))
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (μ/γ^2-(1-μ)/(1-γ)^2-γ-μ+1)/(-2*μ/γ^3-2*(1-μ)/(1-γ)^3-1)
            count += 1
        end
        if count >= maxCount
            Logging.@error "Newton-Raphson did not converge for L1" μ γ iterations=count
            throw(ErrorException("Could not converge on L1 location after $mxCount iterations"))
        end
        q[1] = 1-μ-γ
    elseif point == 2
        Logging.@debug "Iterating on L2 position" μ

        # L2 point beyond secondary
        # Initial guess from Hill's sphere approximation
        γ = cbrt(μ/(3(1-μ)))
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (-μ/γ^2-(1-μ)/(1+γ)^2+γ-μ+1)/(2*μ/γ^3+2*(1-μ)/(1+γ)^3+1)
            count += 1
        end
        if count >= maxCount
            Logging.@error "Newton-Raphson did not converge for L2" μ γ iterations=count
            throw(ErrorException("Could not converge on L2 location after $mxCount iterations"))
        end
        q[1] = 1-μ+γ
    elseif point == 3
        Logging.@debug "Iterating on L3 position" μ

        # L3 point beyond primary
        # Initial guess from first-order series approximation
        γ = 1-7*μ/12
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (μ/(-1-γ)^2+(1-μ)/γ^2-γ-μ)/(-2*μ/(1+γ)^3-2*(1-μ)/γ^3-1)
            count += 1
        end
        if count >= maxCount
            Logging.@error "Newton-Raphson did not converge for L3" μ γ iterations=count
            throw(ErrorException("Could not converge on L3 location after $mxCount iterations"))
        end
        q[1] = -μ-γ
    else
        # L4/5 triangular points, L4 leads secondary, L5 trails it
        # Closed-form solution, each forms equilateral triangle with two primaries
        q[1] = 0.5-μ
        q[2] = (point == 4 ? sin(π/3) : -sin(π/3))
    end

    Logging.@debug "Returning equilibrium point state" point q

    return q
end

"""
    getHamiltonian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64}) -> Float64

Return CR3BP Jacobi constant

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q::AbstractVector{Float64}`: State vector [ndim]

Returns
- `Float64`: Jacobi constant [ndim]

Errors
- Throws `ArgumentError` if state vector is non-finite or too short

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing Jacobi constant,
    or when Jacobi constant is returned

Example
```
H::Float64 = getHamiltonian(dynamicsModel, q)
```
"""
function getHamiltonian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered getHamiltonian (CR3BP)" dynamicsModel q

    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    # Validate state vector length
    n_in::Int64 = length(q)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if n_in < n_simple
        Logging.@error "State vector is too short" required=n_simple actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_simple to calculate Jacobi constant"))
    end

    Logging.@debug "Computing Jacobi constant" q

    U::Float64 = getPseudopotential(dynamicsModel, q)
    v2::Float64 = q[4]^2+q[5]^2+q[6]^2
    JC::Float64 = 2*U-v2

    Logging.@debug "Returning Jacobi constant" JC

    return JC
end

"""
    getJacobiConstant(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64}) -> Float64
    
Return CR3BP Jacobi constant

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q::AbstractVector{Float64}`: State vector [ndim]

Returns
- `Float64`: Jacobi constant [ndim]

Errors
- Throws `ArgumentError` if state vector is too short

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when forwarding to
    `getHamiltonian()`

Notes
- This is a common alternate name for `getHamiltonian` and forwards to that
    method

Example
```
JC::Float64 = getJacobiConstant(dynamicsModel, q)
```
"""
function getJacobiConstant(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered getJacobiConstant (CR3BP)" dynamicsModel q

    # Validate state vector length
    n_in::Int64 = length(q)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if n_in < n_simple
        Logging.@error "State vector is too short" required=n_simple actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_simple to calculate Jacobi constant"))
    end

    Logging.@debug "Forwarding to getHamiltonian" dynamicsModel q[1:n_simple]
    return getHamiltonian(dynamicsModel, q[1:n_simple])
end

"""
    getMassRatios(dynamicsModel::CR3BPDynamicsModel) -> Float64
    
Return CR3BP mass ratio

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object

Returns
- `Float64`: Ratio of secondary to total mass

Errors
- Throws `DomainError` if gravitational parameters are not finite or
    non-positive

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when mass ratio is returned

Example
```
μ::Float64 = getMassRatios(dynamicsModel)
```
"""
function getMassRatios(dynamicsModel::CR3BPDynamicsModel)::Float64
    Logging.@debug "Entered getMassRatios (CR3BP)" dynamicsModel

    # Extract gravitational parameters for both primaries
    μ_1::Float64 = dynamicsModel.primaryData[1].μ
    μ_2::Float64 = dynamicsModel.primaryData[2].μ
    for (label, μ) in (("μ_1", μ_1), ("μ_2", μ_2))
        if !isfinite(μ)
            Logging.@error "Gravitational parameter is not finite" param=label value=μ
            throw(DomainError(μ, "$label must be finite"))
        end
        if μ <= 0.0
            Logging.@error "Gravitational parameter is non-positive" param=label value=μ
            throw(DomainError(μ, "$label must be positive"))
        end
    end

    # Calculate mass ratio
    μ::Float64 = μ_2/(μ_1+μ_2)

    Logging.@debug "Returning mass ratio" massRatio=μ

    return μ
end

"""
    getParameterDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64}) -> Matrix{Float64}

Return derivative of state with respect to parameters

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q_full::Vector{Float64}`: Full state vector [ndim]

Returns
- `Matrix{Float64}`: Derivative of state with respect to parameters [ndim]

Errors
- Throws `ArgumentError` if length of `q_full` is incorrect

Logging
- Emits `@error` logs for thrown errors
- Emits `@info` logs when empty matrix is returned
- Emits `@debug` logs when parameter dependency matrix is returned

Example
```
dqdp::Matrix{Float64} = getParameterDependencies(dynamicsModel, q_full)
```
"""
function getParameterDependencies(dynamicsModel::CR3BPDynamicsModel, q_full::Vector{Float64})::Matrix{Float64}
    Logging.@debug "Entered getParameterDependencies (CR3BP)" dynamicsModel

    # Validate state vector length
    n_full::Int64 = getStateSize(dynamicsModel, FULL)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if length(q_full) != n_full
        Logging.@error "Input state vector has incorrect length" expected=n_full actual=length(q_full)
        throw(ArgumentError("State vector length is $(length(q_full)), but should be $n_full"))
    end

    # CR3BP has no free parameters
    Logging.@info "Returning empty parameter dependency matrix"
    Logging.@debug "Returning parameter dependency matrix" rows=n_simple cols=0
    
    return zeros(Float64, (n_simple,0))
end

"""
    getPrimaryState(dynamicsModel::CR3BPDynamicsModel, primary::Int64) -> Vector{Float64}

Return CR3BP primary state

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `primary::Int64`: Primary index

Returns
- `Vector{Float64}`: Primary state [ndim]

Errors
- Throws `ArgumentError` if `primary` is not between 1 and 2

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing primary state,
    or when primary state is returned

Example
```
q_1::Vector{Float64} = getPrimaryState(dynamicsModel, 1)
```
"""
function getPrimaryState(dynamicsModel::CR3BPDynamicsModel, primary::Int64)::Vector{Float64}
    Logging.@debug "Entered getPrimaryState (CR3BP)" dynamicsModel primary

    # Validate primary index
    if !(1 <= primary <= 2)
        Logging.@error "Invalid primary index" primary
        throw(ArgumentError("Primary must be between 1 and 2, got $primary"))
    end

    μ::Float64 = getMassRatios(dynamicsModel)

    Logging.@debug "Computing primary state" primary μ

    # Both primaries are stationary by definition
    q::Vector{Float64} = zeros(Float64, 6)
    q[1] = (primary == 1 ? -μ : 1-μ)

    Logging.@debug "Returning primary state" primary q

    return q
end

"""
    getPseudopotential(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64}) -> Float64

Return CR3BP pseudo-potential

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q::AbstractVector{Float64}`: State vector [ndim]

Returns
- `Float64`: Pseudo-potential [ndim]

Errors
- Throws `ArgumentError` if state vector is non-finite or too short
- Throws `DomainError` if at location of either primary

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing pseudo-
    potential, or when pseudo-potential is returned

Example
```
U::Float64 = getPseudopotential(dynamicsModel, q)
```
"""
function getPseudopotential(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64})::Float64
    Logging.@debug "Entered getPseudopotential (CR3BP)" dynamicsModel q

    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    # Validate state vector length
    n_in::Int64 = length(q)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if n_in < n_simple
        Logging.@error "State vector is too short" required=n_simple actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_simple to calculate pseudo-potential"))
    end

    μ::Float64 = getMassRatios(dynamicsModel)

    Logging.@debug "Computing pseudo-potential" q μ

    y2::Float64 = q[2]^2
    z2::Float64 = q[3]^2

    # Distance to each primary
    r_13::Float64 = sqrt((q[1]+μ)^2+y2+z2)
    r_23::Float64 = sqrt((q[1]-1+μ)^2+y2+z2)
    if r_13 < 1E-14
        Logging.@error "Location of primary" r_13
        throw(DomainError(r_13, "Distance to primary must be positive"))
    elseif r_23 < 1E-14
        Logging.@error "Location of secondary" r_23
        throw(DomainError(r_23, "Distance to secondary must be positive"))
    end

    U::Float64 = (1-μ)/r_13+μ/r_23+0.5*(q[1]^2+y2)

    Logging.@debug "Returning pseudo-potential" U

    return U
end

"""
    getPseudopotentialHessian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64}) -> Vector{Float64}

Return pseudo-potential Hessian

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q::AbstractVector{Float64}`: State vector [ndim]

Returns
- `Vector{Float64}`: Pseudo-potential Hessian elements [ndim]

Errors
- Throws `ArgumentError` if state vector is non-finite or too short
- Throws `DomainError` if at location of either primary

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing pseudo-
    potential Hessian, or when pseudo-potential Hessian is returned

Example
```
d2Udr2::Float64 = getPseudopotentialHessian(dynamicsModel, q)
```
"""
function getPseudopotentialHessian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64})::Vector{Float64}
    Logging.@debug "Entered getPseudopotentialHessian (CR3BP)" dynamicsModel

    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    # Validate state vector length
    n_in::Int64 = length(q)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if n_in < n_simple
        Logging.@error "State vector is too short" required=n_simple actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_simple to calculate pseudo-potential"))
    end

    μ::Float64 = getMassRatios(dynamicsModel)

    Logging.@debug "Computing pseudo-potential Hessian" q μ

    # x-displacements from eacjh primary
    x_1::Float64 = q[1]+μ
    x_2::Float64 = q[1]-1+μ

    y2::Float64 = q[2]^2
    z2::Float64 = q[3]^2

    # Distance to each primary
    r_13::Float64 = sqrt(x_1^2+y2+z2)
    r_23::Float64 = sqrt(x_2^2+y2+z2)
    if r_13 < 1E-14
        Logging.@error "Location of primary" r_13
        throw(DomainError(r_13, "Distance to primary must be positive"))
    elseif r_23 < 1E-14
        Logging.@error "Location of secondary" r_23
        throw(DomainError(r_23, "Distance to secondary must be positive"))
    end

    r3_13::Float64 = r_13^3
    r3_23::Float64 = r_23^3
    r5_13::Float64 = r3_13*r_13^2
    r5_23::Float64 = r3_23*r_23^2

    # Factor out repeated composite terms
    A_13::Float64 = (1-μ)/r3_13
    A_23::Float64 = μ/r3_23
    B_13::Float64 = 3*(1-μ)/r5_13
    B_23::Float64 = 3*μ/r5_23
    A::Float64 = A_13+A_23
    B::Float64 = B_13+B_23
    Bx::Float64 = B_13*x_1+B_23*x_2

    # Unique elements of pseudo-potential Hessian
    d2Udr2::Vector{Float64} = Vector{Float64}(undef, 6)
    d2Udr2[1] = 1-A+B_13*x_1^2+B_23*x_2^2
    d2Udr2[2] = 1-A+B*y2
    d2Udr2[3] = -A+B*z2
    d2Udr2[4] = Bx*q[2]
    d2Udr2[5] = Bx*q[3]
    d2Udr2[6]=B*q[2]*q[3]

    Logging.@debug "Returning pseudo-potential Hessian" d2Udr2

    return d2Udr2
end

"""
    getPseudopotentialJacobian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64}) -> Vector{Float64}

Return pseudo-potential Jacobian

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q::AbstractVector{Float64}`: State vector [ndim]

Returns
- `Vector{Float64}`: Pseudo-potential Jacobian [ndim]

Errors
- Throws `ArgumentError` if state vector is non-finite or too short
- Throws `DomainError` if at location of either primary

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered, when computing pseudo-
    potential Jacobian, or when pseudo-potential Jacobianis returned

Example
```
dUdr::Float64 = getPseudopotentialJacobian(dynamicsModel, q)
```
"""
function getPseudopotentialJacobian(dynamicsModel::CR3BPDynamicsModel, q::AbstractVector{Float64})::Vector{Float64}
    Logging.@debug "Entered getPseudopotentialJacobian (CR3BP)" dynamicsModel

    # Input validation
    if !all(isfinite, q)
        Logging.@error "Input state vector contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("State vector must contain only finite values"))
    end

    # Validate state vector length
    n_in::Int64 = length(q)
    n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
    if n_in < n_simple
        Logging.@error "State vector is too short" required=n_simple actual=n_in
        throw(ArgumentError("State vector length $n_in is insufficient; need at least $n_simple to calculate pseudo-potential"))
    end

    μ::Float64 = getMassRatios(dynamicsModel)

    Logging.@debug "Computing pseudo-potential Jacobian" q μ

    # x-displacements from eacjh primary
    x_1::Float64 = q[1]+μ
    x_2::Float64 = q[1]-1+μ

    y2::Float64 = q[2]^2
    z2::Float64 = q[3]^2

    # Distance to each primary
    r_13::Float64 = sqrt(x_1^2+y2+z2)
    r_23::Float64 = sqrt(x_2^2+y2+z2)
    if r_13 < 1E-14
        Logging.@error "Location of primary" r_13
        throw(DomainError(r_13, "Distance to primary must be positive"))
    elseif r_23 < 1E-14
        Logging.@error "Location of secondary" r_23
        throw(DomainError(r_23, "Distance to secondary must be positive"))
    end

    r3_13::Float64 = r_13^3
    r3_23::Float64 = r_23^3

    # Factor out repeated composite terms
    A_13::Float64 = (1-μ)/r3_13
    A_23::Float64 = μ/r3_23
    A::Float64 = A_13+A_23

    # Pseudo-potential Jacobian
    dUdr::Vector{Float64} = Vector{Float64}(undef, 3)
    dUdr[1] = q[1]-A_13*x_1-A_23*x_2
    dUdr[2] = q[2]*(1-A)
    dUdr[3] = -q[3]*A

    Logging.@debug "Returning pseudo-potential Jacobian" dUdr

    return dUdr
end

"""
    getStateSize(dynamicsModel::CR3BPDynamicsModel, equationType::EquationType) -> Int64

Return state vector size for CR3BP equations of motion type

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `equationType::EquationType`: Equations of motion type

Returns
- `Int64`: State vector size

Errors
- Throws `ArgumentError` if `equationType` is not found in mapping

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered or when state size is returned

Example
```
n_full::Int64 = getStateSize(dynamicsModel, FULL)
```
"""
function getStateSize(dynamicsModel::CR3BPDynamicsModel, equationType::EquationType)::Int64
    Logging.@debug "Entered getStateSize (CR3BP)" dynamicsModel equationType

    # Guard against enumerated values not present in mapping
    if !haskey(_CR3BP_state_sizes, equationType)
        Logging.@error "Unsupported EquationType for CR3BPDynamicsModel" equationType supported=keys(_CR3BP_state_sizes)
        throw(ArgumentError("Unsupported EquationType: $equationType"))
    end

    n_states::Int64 = _CR3BP_state_sizes[equationType]

    Logging.@debug "Returning state size" equationType stateSize=n_states

    return n_states
end
