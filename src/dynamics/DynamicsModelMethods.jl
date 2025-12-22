"""
DynamicsModel methods

Author: Jonathan Richmond
C: 12/22/25
"""


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
