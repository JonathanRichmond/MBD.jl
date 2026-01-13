"""
DynamicsModel methods

Author: Jonathan Richmond
C: 12/22/25
"""


"""
    appendExtraInitialConditions(mod::AbstractDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)

Append extra initial conditions for a given dynamics model and equation type to
the simple initial state vector.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `q0_simple::Vector{Float64}`: A simple initial state vector.
- `outputEquationType::EquationType`: The output equation formulation to use
  (e.g., `SIMPLE`, `STM`).

Returns
- The initial state vector with extra conditions.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- See the dynamics model specialization for a concrete implementation.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q0_simple = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
q0 = appendExtraInitialConditions(model, q0_simple, STM)
```
"""
function appendExtraInitialConditions(mod::AbstractDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)
    Logging.@debug "appendExtraInitialConditions(Abstract) called" type=typeof(mod)

    Logging.@error "appendExtraInitialConditions not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("appendExtraInitialConditions not implemented for abstract DynamicsModel"))
end

"""
    getCharLengths(mod::AbstractDynamicsModel)

Compute the characteristic length scales for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- The characteristic length scales for the model.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the secondary body's orbital radius.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
lstar = getCharLengths(model)
```
"""
function getCharLengths(mod::AbstractDynamicsModel)
    Logging.@debug "getCharLengths(Abstract) called" type=typeof(mod)

    Logging.@error "getCharLengths not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getCharLengths not implemented for abstract DynamicsModel"))
end

"""
    getCharMasses(mod::AbstractDynamicsModel)

Compute the characteristic mass scales for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- The characteristic mass scales for the model.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the sum of GM values divided by the
  gravitational constant.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
mstar = getCharMasses(model)
```
"""
function getCharMasses(mod::AbstractDynamicsModel)
    Logging.@debug "getCharMasses(Abstract) called" type=typeof(mod)

    Logging.@error "getCharMasses not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getCharMasses not implemented for abstract DynamicsModel"))
end

"""
    getCharTimes(mod::AbstractDynamicsModel)

Compute the characteristic time scales for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- The characteristic time scales for the model.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns sqrt(l*³/ΣGM) where l* is characteristic length
  and ΣGM is the sum of GM values for the primaries.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
tstar = getCharTimes(model)
```
"""
function getCharTimes(mod::AbstractDynamicsModel)
    Logging.@debug "getCharTimes(Abstract) called" type=typeof(mod)

    Logging.@error "getCharTimes not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getCharTimes not implemented for abstract DynamicsModel"))
end

"""
    getMassRatios(mod::AbstractDynamicsModel)

Compute the mass ratios for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- The mass ratios for the model.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the ratio of the secondary body's mass to
  the total mass of the primary bodies.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
μ = getMassRatios(model)
```
"""
function getMassRatios(mod::AbstractDynamicsModel)
    Logging.@debug "getMassRatios(Abstract) called" type=typeof(mod)

    Logging.@error "getMassRatios not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getMassRatios not implemented for abstract DynamicsModel"))
end

"""
    getStateSize(mod::AbstractDynamicsModel, equationType::EquationType)

Return the size of the state vector for a dynamics model and equation type.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `equationType::EquationType`: The equation formulation to use (e.g., `SIMPLE`, `STM`).

Returns
- The number of state variables for the given model and equation type.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- See the dynamics model specialization for a concrete mapping.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
n = getStateSize(model, SIMPLE)
```
"""
function getStateSize(mod::AbstractDynamicsModel, equationType::EquationType)
    Logging.@debug "getStateSize(Abstract) called" type=typeof(mod)

    Logging.@error "getStateSize not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getStateSize not implemented for abstract DynamicsModel"))
end

"""
    shallowClone(mod::AbstractDynamicsModel)

Create a shallow clone of a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance to clone.

Returns
- A new dynamics model instance with copied/reused data as appropriate for the type.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, copies the `primaryData` vector but preserves references
  to individual `BodyData` objects.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
clone = shallowClone(model)
```
"""
function shallowClone(mod::AbstractDynamicsModel)
    Logging.@debug "shallowClone(Abstract) called" type=typeof(mod)

    Logging.@error "shallowClone not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("shallowClone not implemented for abstract DynamicsModel"))
end


"""
    getNumPrimaries(mod::AbstractDynamicsModel) -> Int64

Return the number of primary bodies in a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- `Int64`: The count of primary bodies.

Errors
- Throws `ArgumentError` if `mod` is not a valid `AbstractDynamicsModel`.
- Throws `ErrorException` if `primaryData` is missing or not a vector.
- Throws `ArgumentError` if model-specific constraints are violated (e.g., CR3BP
  requires exactly 2 primaries).

Notes
- Validates that the model instance is well-formed.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
n = getNumPrimaries(model)
```
"""
function getNumPrimaries(mod::AbstractDynamicsModel)
    Logging.@debug "getNumPrimaries(Abstract) called" type=typeof(mod)

    if mod === nothing || !isa(mod, AbstractDynamicsModel)
        Logging.@error "Invalid model supplied to getNumPrimaries" type=typeof(mod)
        throw(ArgumentError("mod must be an AbstractDynamicsModel instance"))
    end
    if !hasproperty(mod, :primaryData)
        Logging.@error "Model has no primaryData field" type=typeof(mod)
        throw(ErrorException("model does not define primaryData"))
    end

    prim = getfield(mod, :primaryData)
    if !(isa(prim, AbstractVector))
        Logging.@error "primaryData is not a vector" typeof_primaryData=typeof(prim)
        throw(ArgumentError("primaryData must be a vector"))
    end

    n::Int64 = length(prim)

    Logging.@info "Computed number of primaries" count=n
    return n
end


"""
    appendExtraInitialConditions(mod::CR3BPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType) -> Vector{Float64}

Append extra initial conditions for a given equation type to a simple initial state vector.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `q0_simple::Vector{Float64}`: A simple initial state vector.
- `outputEquationType::EquationType`: The output equation formulation (e.g., `SIMPLE`, `STM`, `ARCLENGTH`).

Returns
- `Vector{Float64}`: The initial state vector expanded to the required size for the output equation type.
  - For `SIMPLE`: returns the input vector unchanged.
  - For other types: returns the input vector followed by additional initial conditions.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `q0_simple` is not a non-empty `Vector{Float64}`.
- Throws `ArgumentError` if the input vector size does not match the SIMPLE state size.

Notes
- Supported output equation types: `SIMPLE`, `STM`, `ARCLENGTH`, `MOMENTUM`, `FULL`.
- When expanding to STM or higher, diagonal identity-like patterns may be used for
  state transition matrix initialization.
- The function emits `Logging` messages at `@debug`, `@error`, and `@info` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q0_simple = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
q0 = appendExtraInitialConditions(model, q0_simple, STM)
```
"""
function appendExtraInitialConditions(mod::CR3BPDynamicsModel, q0_simple::Vector{Float64}, outputEquationType::EquationType)
    Logging.@debug "appendExtraInitialConditions(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to appendExtraInitialConditions" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !isa(q0_simple, AbstractVector)
        Logging.@error "q0_simple must be a vector" type=typeof(q0_simple)
        throw(ArgumentError("q0_simple must be a vector"))
    end
    if isempty(q0_simple)
        Logging.@error "q0_simple is empty"
        throw(ArgumentError("q0_simple must be non-empty"))
    end

    n_in::Int16 = Int16(length(q0_simple))
    n_simple::Int16 = Int16(getStateSize(mod, SIMPLE))
    if n_in != n_simple
        Logging.@error "Input state vector size does not match SIMPLE size" n_in=n_in expected=n_simple
        throw(ArgumentError("Input state vector size ($(n_in)) does not match SIMPLE state size ($(n_simple))"))
    end

    n_out::Int16 = Int16(getStateSize(mod, outputEquationType))
    q0::Vector{Float64} = zeros(Float64, n_out)
    q0[1:n_simple] = q0_simple
    if n_out > n_simple
        n_STM::Int16 = Int16(getStateSize(mod, STM))
        [q0[j] = 1.0 for j in n_simple+1:n_simple+1:n_STM]
    end

    Logging.@info "Appended extra initial conditions" n_in=n_in n_out=n_out
    return q0
end

"""
    getCharLengths(mod::CR3BPDynamicsModel) -> Float64

Compute the characteristic length scale for a CR3BP dynamics model.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `Float64`: characteristic length, in km.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if the secondary's orbital radius is not finite or positive.

Notes
- Characteristic length is the semi-major axis (circular radius) of the secondary.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
lstar = getCharLengths(model)
```
"""
function getCharLengths(mod::CR3BPDynamicsModel)
    Logging.@debug "getCharLengths(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getCharLengths" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    lstar::Float64 = mod.primaryData[2].a
    if !(isfinite(lstar) && (lstar > 0))
        Logging.@error "Invalid characteristic length (a) for secondary" a=lstar
        throw(ArgumentError("Characteristic length must be finite and positive"))
    end

    Logging.@info "Computed characteristic length" lstar=lstar
    return lstar
end

"""
    getCharMasses(mod::CR3BPDynamicsModel) -> Float64

Compute the characteristic mass for a CR3BP dynamics model.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `Float64`: The characteristic mass, in kg.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if the computed mass is not finite or positive.

Notes
- Characteristic mass is the sum of gravitational parameters divided by the
  gravitational constant.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
mstar = getCharMasses(model)
```
"""
function getCharMasses(mod::CR3BPDynamicsModel)
    Logging.@debug "getCharMasses(CR3BP) called"
    
    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getCharMasses" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    mstar::Float64 = (mod.primaryData[1].μ + mod.primaryData[2].μ) / GRAVITY
    if !(isfinite(mstar) && (mstar > 0))
        Logging.@error "Invalid characteristic mass" mass=mstar
        throw(ArgumentError("Characteristic mass must be finite and positive"))
    end

    Logging.@info "Computed characteristic mass" mass=mstar
    return mstar
end

"""
    getCharTimes(mod::CR3BPDynamicsModel) -> Float64

Compute the characteristic time scale for a CR3BP dynamics model.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `Float64`: The characteristic time, in seconds.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `l*` or ΣGM are not finite/positive, or if tstar is invalid.

Notes
- Characteristic time is derived from the characteristic length and total GM.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
tstar = getCharTimes(model)
```
"""
function getCharTimes(mod::CR3BPDynamicsModel)
    Logging.@debug "getCharTimes(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getCharTimes" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    lstar::Float64 = getCharLengths(mod)
    totalμ::Float64 = mod.primaryData[1].μ + mod.primaryData[2].μ
    if !(isfinite(totalμ) && (totalμ > 0))
        Logging.@error "Invalid total GM for primaries" totalμ=totalμ
        throw(ArgumentError("Total GM must be finite and positive"))
    end
    tstar::Float64 = sqrt(lstar^3 / totalμ)
    if !(isfinite(tstar) && (tstar > 0))
        Logging.@error "Invalid characteristic time" tstar=tstar
        throw(ArgumentError("Characteristic time must be finite and positive"))
    end

    Logging.@info "Computed characteristic time" tstar=tstar
    return tstar
end

"""
    getMassRatios(mod::CR3BPDynamicsModel) -> Float64
    
Compute the mass ratio for a CR3BP dynamics model.
    
Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `Float64`: The mass ratio, where 0 < μ < 1.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if computed mass ratio is invalid (not finite, outside (0,1)).
  
Notes
- Mass ratio represents the secondary body's fractional contribution to total mass.
- Used in CR3BP analysis for nondimensionalization and stability calculations.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.
  
Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
μ = getMassRatios(model)
```
"""
function getMassRatios(mod::CR3BPDynamicsModel)
    Logging.@debug "getMassRatios(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getMassRatios" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    totalμ::Float64 = mod.primaryData[1].μ + mod.primaryData[2].μ
    if !(isfinite(totalμ) && (totalμ > 0))
        Logging.@error "Invalid GM values for mass ratio" μ1=mod.primaryData[1].μ μ2=mod.primaryData[2].μ totalμ=totalμ
        throw(ArgumentError("GM values must be finite and positive"))
    end
    μ::Float64 = mod.primaryData[2].μ / totalμ
    if !(isfinite(μ) && (μ > 0) && (μ < 1))
        Logging.@error "Invalid mass ratio computed" μ=μ
        throw(ArgumentError("Mass ratio must be finite and positive"))
    end

    Logging.@info "Computed mass ratio" μ=μ
    return μ
end

"""
    getStateSize(mod::CR3BPDynamicsModel, equationType::EquationType) -> Int64

Return the size of the state vector for a CR3BP dynamics model and equation type.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `equationType::EquationType`: The equation formulation to use (e.g., `SIMPLE`, `STM`).

Returns
- `Int64`: The number of state variables for a CR3BP dynamics model and given equation type.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
    contain exactly 2 primary bodies.
- Throws `ArgumentError` if `equationType` is unsupported for CR3BP.
- Throws `ArgumentError` if the computed state size is not a positive integer.

Notes
- Supported `EquationType → size` mapping:
    - `SIMPLE → 6`
    - `STM → 42`
    - `ARCLENGTH → 43`
    - `MOMENTUM → 43`
    - `FULL → 44`
- Returned value is normalized to `Int64` and validated.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
n = getStateSize(model, SIMPLE)
```
"""
function getStateSize(mod::CR3BPDynamicsModel, equationType::EquationType)
    Logging.@debug "getStateSize(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getStateSize" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    sizeMap = Dict(SIMPLE => 6, STM => 42, ARCLENGTH => 43, MOMENTUM => 43, FULL => 44)
    if !haskey(sizeMap, equationType)
        Logging.@error "Unsupported equation type for CR3BP" equationType=equationType
        throw(ArgumentError("Unsupported equation type for CR3BPDynamicsModel"))
    end

    size = sizeMap[equationType]
    if !(isa(size, Integer) && (size > 0))
        Logging.@error "Invalid state size computed" equationType=equationType size=size
        throw(ArgumentError("State size must be a positive integer"))
    end
    size64::Int64 = Int64(size)

    Logging.@info "Computed state size for CR3BPDynamicsModel" equationType=equationType size=size64
    return size64
end

"""
    shallowClone(mod::CR3BPDynamicsModel) -> CR3BPDynamicsModel

Create a shallow clone of a `CR3BPDynamicsModel` instance.

Arguments
- `mod::CR3BPDynamicsModel`: Valid `CR3BPDynamicsModel` to clone.

Returns
- `CR3BPDynamicsModel`: A new `CR3BPDynamicsModel` with a copied `primaryData` vector.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel`.

Behavior
- - Copies `primaryData` vector.

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
clone = shallowClone(model)
```
"""
function shallowClone(mod::CR3BPDynamicsModel)
    Logging.@debug "shallowClone(CR3BP) called"

    if mod === nothing || !isa(mod, MBD.CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to shallowClone" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end

    clone = CR3BPDynamicsModel(copy(mod.primaryData))

    Logging.@info "Created shallow clone of CR3BPDynamicsModel" n_bodies=length(clone.primaryData)
    return clone
end
