"""
DynamicsModel methods

Author: Jonathan LeFevre Richmond
C: 5/4/26
U: 5/6/26
"""


"""
    adjustInitialConditions(dynamicsModel::AbstractDynamicsModel, q0::Vector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)

Return initial conditions for output equations of motion type

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q0::Vector{Float64}`: Initial conditions
- `inputEquationType::EquationType`: Equations of motion type for `q0`
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector::Float64`: Initial conditions

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
function adjustInitialConditions(dynamicsModel::AbstractDynamicsModel, q0::Vector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered generic adjustInitialConditions" dynamicsModel inputEquationType outputEquationType

    Logging.@error "adjustInitialConditions is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(adjustInitialConditions, (dynamicsModel, q0, inputEquationType, outputEquationType)))
end

"""
    appendExtraInitialConditions(dynamicsModel::AbstractDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)

Return initial conditions for output equations of motion type

Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
- `q0_simple::Vector{Float64}`: Simple initial conditions
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
- This is a common special case of `adjustInitialConditions` and forwards to
    that method
- For `CR3BPDynamicsModel`, the only non-zero values are in the simple state
    and STM

Example
```
q0_full::Vector{Float64} = appendExtraInitialConditions(dynamicsModel, q0_simple, FULL)
```
"""
function appendExtraInitialConditions(dynamicsModel::AbstractDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered generic appendExtraInitialConditions" dynamicsModel outputEquationType

    Logging.@error "appendExtraInitialConditions is not implemented for this dynamics model type" type=typeof(dynamicsModel)
    throw(MethodError(appendExtraInitialConditions, (dynamicsModel, q0_simple, outputEquationType)))
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
    getStateSize(dynamicsModel::AbstractDynamicsModel, equationType::EquationType)

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
    adjustInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0::Vector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)

Return initial conditions for CR3BP output equations of motion type

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q0::Vector{Float64}`: Initial conditions [ndim]
- `inputEquationType::EquationType`: Equations of motion type for `q0`
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector::Float64`: Initial conditions [ndim]

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
function adjustInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0::Vector{Float64}, inputEquationType::EquationType, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered appendExtraInitialConditions" dynamicsModel inputEquationType outputEquationType

    # Validate input state vector length against expected size for inputEquationType
    n_in::Int64 = getStateSize(dynamicsModel, inputEquationType)
    if length(q0) != n_in
        Logging.@error "Input state vector has incorrect length" expected=n_in actual=length(q0)
        throw(ArgumentError("State vector length is $(length(q0)), but should be $n_in"))
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

        # Intiialize STM block as flattened identity matrix only when input doesn't already contain STM block but output does
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
    appendExtraInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)

Return initial conditions for CR3BP output equations of motion type

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `q0_simple::Vector{Float64}`: Simple initial conditions [ndim]
- `outputEquationType::EquationType`: Output equations of motion type

Returns
- `Vector{Float64}`: Initial conditions [ndim]

Errors
- No additional error checking

Logging
- Emits `@debug` logs when function is entered or when forwarding to
    `adjustInitialConditions`

Notes
- This is a common special case of `adjustInitialConditions` and forwards to
    that method

Example
```
q0_full::Vector{Float64} = appendExtraInitialConditions(dynamicsModel, q0_simple, FULL)
```
"""
function appendExtraInitialConditions(dynamicsModel::CR3BPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)::Vector{Float64}
    Logging.@debug "Entered appendExtraInitialConditions" dynamicsModel outputEquationType

    Logging.@debug "Forwarding to adjustInitialConditions" dynamicsModel outputEquationType
    return adjustInitialConditions(dynamicsModel, q0_simple, SIMPLE, outputEquationType)
end

"""
    getCharLengths(dynamicsModel::CR3BPDynamicsModel)
    
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
    getCharMasses(dynamicsModel::CR3BPDynamicsModel)
    
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

    # Extract gravitational parameterts for both primaries
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
    getCharTimes(dynamicsModel::CR3BPDynamicsModel)
    
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

    # Extract characteristic length and gravitational parameter
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
    getStateSize(dynamicsModel::CR3BPDynamicsModel, equationType::EquationType)

Return state vector size for CR3BP equations of motion type

Arguments
- `dynamicsModel::CR3BPDynamicsModel`: `CR3BPDynamicsModel` object
- `equationType::EquationType`: Equations of motion type

Returns
- `Int64`: State vector size

Errors
- Throws `ArgumentError` if equationType is not found in mapping

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
