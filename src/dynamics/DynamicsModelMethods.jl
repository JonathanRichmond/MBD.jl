"""
DynamicsModel methods

Author: Jonathan LeFevre Richmond
C: 5/4/26
U: 5/5/26
"""


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
