"""
DynamicsModel methods

Author: Jonathan LeFevre Richmond
C: 12/22/25
U: 1/15/26

ON HOLD: checkSTM (Propagator), evaluateEquations (computeDerivatives)

TO DO: getTidalAcceleration, get2BApproximation, primaryInertial2Rotating,
    rotating2PrimaryEclipJ2000, rotating2PrimaryInertial,
    rotating2SunEclipJ2000, secondaryEclipJ20002Rotating
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
    extractStateTransitionMatrix(mod::AbstractDynamicsModel, q::Vector{Float64}) -> Matrix{Float64}

Extract the state transition matrix from an STM state vector.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `q::Vector{Float64}`: An STM state vector containing the state transition matrix block in column-major order.

Returns
- `Matrix{Float64}`: The state transition matrix extracted from the state vector.

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
q = appendExtraInitialConditions(model, q_simple, STM)  # State with STM
Φ = extractStateTransitionMatrix(model, q)
```
"""
function extractStateTransitionMatrix(mod::AbstractDynamicsModel, q::Vector{Float64})
    Logging.@debug "extractStateTransitionMatrix(Abstract) called" type=typeof(mod)

    Logging.@error "extractStateTransitionMatrix not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("extractStateTransitionMatrix not implemented for abstract DynamicsModel"))
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
    getDistance2Primary(mod::AbstractDynamicsModel, primaryID::Int64, pos::Vector{Float64}) -> Float64

Compute the distance from a position to a primary body in a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `primaryID::Int64`: The ID of the primary body.
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Float64`: The Euclidean distance from the position to the primary body.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
d = getDistance2Primary(model, 1, [1.0, 0, 0])
```
"""
function getDistance2Primary(mod::AbstractDynamicsModel, primaryID::Int64, pos::Vector{Float64})
    Logging.@debug "getDistance2Primary(Abstract) called" type=typeof(mod)

    Logging.@error "getDistance2Primary not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getDistance2Primary not implemented for abstract DynamicsModel"))
end

"""
    getEnergy(mod::AbstractDynamicsModel, q::Vector{Float64})

Compute the energy of a state in a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `q::Vector{Float64}`: A 6D state vector in normalized coordinates.

Returns
- `Float64`: The energy of the state.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the Jacobi constant value for the state.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
E = getEnergy(model, [1.0, 0, 0, 0, 0.1, 0])
```
"""
function getEnergy(mod::AbstractDynamicsModel, q::Vector{Float64})
    Logging.@debug "getEnergy(Abstract) called" type=typeof(mod)

    Logging.@error "getEnergy not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getEnergy not implemented for abstract DynamicsModel"))
end

"""
    getEquationsOfMotion(mod::AbstractDynamicsModel) -> AbstractEquationsOfMotion

Construct the equations of motion for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- `AbstractEquationsOfMotion`: An equations of motion instance encapsulating the dynamics equations for the model.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
eom = getEquationsOfMotion(model)
```
"""
function getEquationsOfMotion(mod::AbstractDynamicsModel)
    Logging.@debug "getEquationsOfMotion(Abstract) called" type=typeof(mod)

    Logging.@error "getEquationsOfMotion not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getEquationsOfMotion not implemented for abstract DynamicsModel"))
end

"""
    getEquilibriumPoint(mod::AbstractDynamicsModel, pointID::Int64)

Compute the position of an equilibrium point in a dynamics model rotating frame.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `pointID::Int64`: The ID of the equilibrium point.

Returns
- `Vector{Float64}`: The 3D position vector of the equilibrium point in normalized coordinates in the rotating frame.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the Lagrange equilibrium point position in
  the CR3BP rotating frame.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
E1 = getEquilibriumPoint(model, 1)
```
"""
function getEquilibriumPoint(mod::AbstractDynamicsModel, pointID::Int64)
    Logging.@debug "getEquilibriumPoint(Abstract) called" type=typeof(mod)

    Logging.@error "getEquilibriumPoint not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getEquilibriumPoint not implemented for abstract DynamicsModel"))
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
    getPrimaryState(mod::AbstractDynamicsModel, primaryID::Int64) -> Vector{Float64}

Compute the state vector of a primary body in a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `primaryID::Int64`: The ID of the primary body.

Returns
- `Vector{Float64}`: The state vector of the primary body in normalized coordinates in the rotating frame.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- For `CR3BPDynamicsModel`, returns the state of primary 1 or 2 in the rotating frame.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q_1 = getPrimaryState(model, 1)
```
"""
function getPrimaryState(mod::AbstractDynamicsModel, primaryID::Int64)
    Logging.@debug "getPrimaryState(Abstract) called" type=typeof(mod)

    Logging.@error "getPrimaryState not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getPrimaryState not implemented for abstract DynamicsModel"))
end

"""
    getPseudopotential(mod::AbstractDynamicsModel, pos::Vector{Float64})

Compute the pseudopotential at a given position in a dynamics model in the rotating frame.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Float64`: The pseudopotential value at the given position.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
U = getPseudopotential(model, [1.0, 0, 0])
```
"""
function getPseudopotential(mod::AbstractDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotential(Abstract) called" type=typeof(mod)

    Logging.@error "getPseudopotential not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getPseudopotential not implemented for abstract DynamicsModel"))
end

"""
    getPseudopotentialGradient(mod::AbstractDynamicsModel, pos::Vector{Float64})

Compute the gradient of the pseudopotential at a given position in a dynamics model in the rotating frame.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Vector{Float64}`: The derivatives of the pseudopotential value at the given position.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
dU = getPseudopotentialGradient(model, [1.0, 0, 0])
```
"""
function getPseudopotentialGradient(mod::AbstractDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotentialGradient(Abstract) called" type=typeof(mod)

    Logging.@error "getPseudopotentialGradient not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getPseudopotentialGradient not implemented for abstract DynamicsModel"))
end

"""
    getPseudopotentialHessian(mod::AbstractDynamicsModel, pos::Vector{Float64})

Compute the Hessian of the pseudopotential at a given position in a dynamics model in the rotating frame.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Vector{Float64}`: The Jacobian of the pseudopotential value at the given position.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
ddU = getPseudopotentialHessian(model, [1.0, 0, 0])
```
"""
function getPseudopotentialHessian(mod::AbstractDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotentialHessian(Abstract) called" type=typeof(mod)
    
    Logging.@error "getPseudopotentialHessian not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("getPseudopotentialHessian not implemented for abstract DynamicsModel"))
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
    isEpochIndependent(mod::AbstractDynamicsModel) -> Bool

Determine whether a dynamics model's equations of motion are epoch-independent.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.

Returns
- `Bool`: True if the model is epoch-independent (time-invariant); false otherwise.

Errors
- Throws `ErrorException` if not implemented for the model type.

Notes
- This is a generic method that dispatches on model type.
- Epoch-independent models have equations of motion that do not explicitly depend on time.
- The function emits `Logging` messages at `@debug` and `@error` levels to aid
  troubleshooting.
- See the dynamics model specialization for a concrete implementation.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
indep = isEpochIndependent(model)
```
"""
function isEpochIndependent(mod::AbstractDynamicsModel)
    Logging.@debug "isEpochIndependent(Abstract) called" type=typeof(mod)

    Logging.@error "isEpochIndependent not implemented for abstract DynamicsModel" type=typeof(mod)
    throw(ErrorException("isEpochIndependent not implemented for abstract DynamicsModel"))
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
    getEpochDependencies(mod::AbstractDynamicsModel, q_full::Vector{Float64}) -> Matrix{Float64}

Extract epoch sensitivities from a full state vector for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `q_full::Vector{Float64}`: A full state vector including epoch dependencies.

Returns
- `Matrix{Float64}`: A matrix of shape (n_simple, 1) containing epoch sensitivities.
  For epoch-independent models, returns an (n_simple, 0) empty matrix.
  For epoch-dependent models, returns sensitivities ∂q/∂epoch.

Errors
- Throws `ArgumentError` if `mod` is not a valid `AbstractDynamicsModel`.
- Throws `ArgumentError` if `q_full` is not a `Vector{Float64}` or is empty.
- Throws `ArgumentError` if `q_full` size does not match `getStateSize(mod, FULL)`.
- Throws `ArgumentError` if `q_full` contains non-finite values.
- Throws `ErrorException` if computed matrix contains non-finite values.

Notes
- State vector structure: [simple | STM | epoch_deps (variable)]
- Epoch dependencies start after STM block at index n_STM+1.
- For epoch-independent models (like CR3BP), returns empty (n_simple, 0) matrix.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q_full = appendExtraInitialConditions(model, [1.0, 0, 0, 0, 0.1, 0], FULL)
∂q∂E = getEpochDependencies(model, q_full)
```
"""
function getEpochDependencies(mod::AbstractDynamicsModel, q_full::Vector{Float64})
    Logging.@debug "getEpochDependencies(Abstract) called"

    if mod === nothing || !isa(mod, AbstractDynamicsModel)
        Logging.@error "Invalid model supplied to getEpochDependencies" type=typeof(mod)
        throw(ArgumentError("mod must be an AbstractDynamicsModel instance"))
    end
    if !(isa(q_full, AbstractVector) && all(isa.(q_full, Float64)))
        Logging.@error "q_full must be a Vector{Float64}" type=typeof(q_full)
        throw(ArgumentError("q_full must be a Vector{Float64}"))
    end
    if isempty(q_full)
        Logging.@error "q_full is empty"
        throw(ArgumentError("q_full must be non-empty"))
    end
    if !all(isfinite, q_full)
        non_finite_count = count(!isfinite, q_full)
        Logging.@error "q_full contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("q_full must contain only finite values"))
    end

    n_in::Int16 = Int16(length(q_full))
    n_full::Int16 = Int16(getStateSize(mod, FULL))
    if n_in != n_full
        Logging.@error "Input state vector size does not match FULL size" n_in=n_in expected=n_full
        throw(ArgumentError("Input state vector size ($(n_in)) does not match FULL state size ($(n_full))"))
    end
    n_STM::Int16 = Int16(getStateSize(mod, STM))
    n_simple::Int16 = Int16(getStateSize(mod, SIMPLE))

    # Return empty matrix if no epoch dependencies (epoch-independent model)
    if n_STM == n_full
        Logging.@info "No epoch dependencies present (epoch-independent model)" result_size=(n_simple, 0)
        return zeros(Float64, (n_simple,0))
    end

    Logging.@debug "Extracting epoch dependencies" n_simple=n_simple epoch_block_start=(n_STM+1) epoch_block_end=n_full

    # Extract epoch dependencies from q_full
    ∂q∂E::Matrix{Float64} = reshape(q_full[n_STM+1:end], (n_simple,1))
    Logging.@debug "Epoch dependency matrix extracted" shape=size(∂q∂E)
    if !all(isfinite, ∂q∂E)
        non_finite_count = count(!isfinite, ∂q∂E)
        Logging.@error "Extracted epoch dependencies contain non-finite values" non_finite_count=non_finite_count
        throw(ErrorException("Extracted epoch dependency matrix contains non-finite values"))
    end

    Logging.@info "Successfully extracted epoch dependencies" size=size(∂q∂E)
    return ∂q∂E
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
    getParameterDependencies(mod::AbstractDynamicsModel, q_full::Vector{Float64}) -> Matrix{Float64}

Extract parameter sensitivities from a full state vector for a dynamics model.

Arguments
- `mod::AbstractDynamicsModel`: A dynamics model instance.
- `q_full::Vector{Float64}`: A full state vector including parameter dependencies.

Returns
- `Matrix{Float64}`: A matrix of shape (n_simple, n_params) where:
  - Rows correspond to simple state elements
  - Columns correspond to parameter dependencies
  - Each column [r,c] is the derivative of state r with respect to parameter c.
  - Returns empty matrix if no parameter dependencies are present.

Errors
- Throws `ArgumentError` if `mod` is not a valid `AbstractDynamicsModel`.
- Throws `ArgumentError` if `q_full` is not a `Vector{Float64}` or is empty.
- Throws `ArgumentError` if `q_full` size does not match `getStateSize(mod, FULL)`.
- Throws `ErrorException` if computed matrix contains non-finite values.

Notes
- State vector structure: [simple | STM | param_deps (variable)]
- Parameter dependencies start after STM block at index n_STM+1.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q_full = [state, STM, dq/dparam]
dqdparam = getParameterDependencies(model, q_full)
```
"""
function getParameterDependencies(mod::CR3BPDynamicsModel, q_full::Vector{Float64})
    Logging.@debug "getParameterDependencies(Abstract) called"

    if mod === nothing || !isa(mod, AbstractDynamicsModel)
        Logging.@error "Invalid model supplied to getParameterDependencies" type=typeof(mod)
        throw(ArgumentError("mod must be an AbstractDynamicsModel instance"))
    end
    if !(isa(q_full, AbstractVector) && all(isa.(q_full, Float64)))
        Logging.@error "q_full must be a Vector{Float64}" type=typeof(q_full)
        throw(ArgumentError("q_full must be a Vector{Float64}"))
    end
    if isempty(q_full)
        Logging.@error "q_full is empty"
        throw(ArgumentError("q_full must be non-empty"))
    end
    if !all(isfinite, q_full)
        non_finite_count = count(!isfinite, q_full)
        Logging.@error "q_full contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("q_full must contain only finite values"))
    end

    n_in::Int16 = Int16(length(q_full))
    n_full::Int16 = Int16(getStateSize(mod, FULL))
    if n_in != n_full
        Logging.@error "Input state vector size does not match FULL size" n_in=n_in expected=n_full
        throw(ArgumentError("Input state vector size ($(n_in)) does not match FULL state size ($(n_full))"))
    end
    n_STM::Int16 = Int16(getStateSize(mod, STM))
    n::Int16 = n_full-n_STM
    n_simple::Int16 = Int16(getStateSize(mod, SIMPLE))
    
    # Return empty matrix if no parameter dependencies
    if n == 0
        Logging.@info "No parameter dependencies present in state vector"
        return zeros(Float64, (n_simple,0))
    end

    # Compute number of parameter sets
    if mod(n, n_simple) != 0
        Logging.@error "Parameter dependencies block size not divisible by n_simple" n=n n_simple=n_simple
        throw(ArgumentError("Parameter dependencies size is inconsistent with state dimension"))
    end
    n_params::Int16 = n/n_simple

    Logging.@debug "Extracting parameter dependencies" n_params=n_params n_simple=n_simple params_block_start=(n_STM+1) params_block_end=n_full
    
    # Extract parameter dependencies from q_full
    ∂q∂param::Matrix{Float64} = zeros(Float64, (n_simple,n_params))
    for r::Int16 in 1:n_simple
        for c::Int16 in 1:n_params
            idx::Int16 = n_STM+n_simple*(c-1)+r
            if (idx < 1) || (idx > n_full)
                Logging.@error "Index out of bounds during extraction" idx=idx length=length(q_full)
                throw(BoundsError("Parameter dependency index out of bounds"))
            end
            ∂q∂param[r,c] = q_full[idx]
            if !isfinite(∂q∂param[r,c])
                Logging.@error "Non-finite parameter dependency extracted" row=r col=c value=∂q∂param[r,c]
                throw(ErrorException("Parameter dependency matrix contains non-finite values"))
            end
        end
    end

    Logging.@info "Computed parameter dependencies matrix" size=(n_simple, n_params) n_params=n_params
    return ∂q∂param
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
    if !(isa(q0_simple, AbstractVector) && all(isa.(q0_simple, Float64)))
        Logging.@error "q0_simple must be a Vector{Float64}" type=typeof(q0_simple)
        throw(ArgumentError("q0_simple must be a Vector{Float64}"))
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
    extractStateTransitionMatrix(mod::CR3BPDynamicsModel, q::Vector{Float64}) -> Matrix{Float64}

Extract the 6x6 state transition matrix from a CR3BP state vector.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `q::Vector{Float64}`: A state vector of at least size 42 (SIMPLE + STM blocks) in column-major order.

Returns
- `Matrix{Float64}`: The 6x6 state transition matrix Φ extracted from indices 7-42.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `q` is not a `Vector{Float64}`, is empty, or is too short
  to contain the STM block.
- Throws `ErrorException` if extracted matrix contains non-finite values.

Notes
- State vector structure: [simple (6) | STM (36)] for minimum size 42 in
  column-major order..
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.
- Useful for extracting linearized dynamics information from propagated states.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q = appendExtraInitialConditions(model, q_simple, STM)  # State with STM
Φ = extractStateTransitionMatrix(model, q)
```
"""
function extractStateTransitionMatrix(mod::CR3BPDynamicsModel, q::Vector{Float64})
    Logging.@debug "extractStateTransitionMatrix(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to extractStateTransitionMatrix" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(q, AbstractVector) && all(isa.(q, Float64)))
        Logging.@error "q must be a Vector{Float64}" type=typeof(q)
        throw(ArgumentError("q must be a Vector{Float64}"))
    end
    if isempty(q)
        Logging.@error "q is empty"
        throw(ArgumentError("q must be non-empty"))
    end
    if !all(isfinite, q)
        Logging.@error "q contains non-finite values" non_finite_count=count(!isfinite, q)
        throw(ArgumentError("q must contain only finite values"))
    end

    n_in::Int16 = Int16(length(q))
    n_STM::Int16 = Int16(getStateSize(mod, STM))
    if n_in < n_STM
        Logging.@error "State vector q is too short to contain STM" length=n_in required=n_STM
        throw(ArgumentError("q must have at least $(n_STM) elements for STM extraction"))
    end
    n_simple::Int16 = Int16(getStateSize(mod, SIMPLE))

    try
        Φ::Matrix{Float64} = reshape(q[n_simple+1:n_STM], (6,6))
        Logging.@debug "State transition matrix extracted" size=(6,6) stm_indices="$(n_simple+1):$(n_STM)"
        if !all(isfinite, Φ)
            non_finite_count = count(!isfinite, Φ)
            Logging.@error "Extracted STM contains non-finite values" non_finite_count=non_finite_count
            throw(ErrorException("Extracted state transition matrix contains non-finite values"))
        end
        frobeniusNorm = LinearAlgebra.norm(Φ)

        Logging.@info "Successfully extracted state transition matrix" size=(6,6) norm_frobenius=frobeniusNorm
        return Φ
    catch e
        if isa(e, ArgumentError) || isa(e, ErrorException)
            Logging.@error "Failed to extract state transition matrix" error=e.msg
            rethrow()
        else
            Logging.@error "Unexpected error during STM extraction" error=string(e)
            throw(ErrorException("Failed to extract state transition matrix: $(e)"))
        end
    end
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
    getDistance2Primary(mod::CR3BPDynamicsModel, primaryID::Int64, pos::Vector{Float64}) -> Float64

Compute the distance from a state to a primary body in the CR3BP rotating frame.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `primaryID::Int64`: The primary body identifier (1 for primary, 2 for secondary).
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Float64`: The Euclidean distance from the position to the primary body.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pos` is not a `Vector{Float64}` or has fewer than 3 elements.
- Throws `ArgumentError` if position contains non-finite values.
- Throws `ArgumentError` if `primaryID` is not 1 or 2.
- Throws `ErrorException` if computed distance is invalid (non-finite or negative).

Notes
- Velocity components (elements 4-6) are ignored for distance calculation.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
d = getDistance2Primary(model, 1, [1.0, 0, 0])
```
"""
function getDistance2Primary(mod::CR3BPDynamicsModel, primaryID::Int64, pos::Vector{Float64})
    Logging.@debug "getDistance2Primary(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getDistance2Primary" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pos, AbstractVector) && all(isa.(pos, Float64)))
        Logging.@error "Position vector pos must be a Vector{Float64}" type=typeof(pos)
        throw(ArgumentError("pos must be a Vector{Float64}"))
    end
    if length(pos) < 3
        Logging.@error "Position vector pos is too short" length=length(pos) required=3
        throw(ArgumentError("pos must have at least 3 elements (position components)"))
    end
    if !all(isfinite, pos)
        non_finite_count = count(!isfinite, pos)
        Logging.@error "pos contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("pos must contain only finite values"))
    end
    if !(isa(primaryID, Integer) && (primaryID == 1 || primaryID == 2))
        Logging.@error "Invalid primaryID supplied" primaryID=primaryID
        throw(ArgumentError("primaryID must be 1 or 2"))
    end

    qP::Vector{Float64} = getPrimaryState(mod, primaryID)

    d::Float64 = sqrt((pos[1]-qP[1])^2+(pos[2]-qP[2])^2+(pos[3]-qP[3])^2)
    if !isfinite(d) || d < 0
        Logging.@error "Invalid distance computed" primaryID=primaryID distance=d
        throw(ErrorException("Distance to primary is non-finite or negative"))
    end

    Logging.@info "Computed distance to primary" primaryID=primaryID distance=d state_position=pos[1:3] primary_position=qP[1:3]
    return d
end

"""
    getEnergy(mod::CR3BPDynamicsModel, q::Vector{Float64}) -> Float64

Compute the Jacobi constant (energy) for a state in a CR3BP dynamics model.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `q::Vector{Float64}`: A 6D state vector in normalized coordinates.

Returns
- `Float64`: The Jacobi constant (CR3BP energy) at the given state.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `q` is not a `Vector{Float64}` or has fewer than 6 elements.
- Throws `ErrorException` if computed energy is non-finite or invalid.

Notes
- The Jacobi constant is a conserved quantity in the CR3BP.
- Lower JC values indicate higher energy and more freedom to traverse the system.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
JC = getEnergy(model, [1.0, 0, 0, 0, 0.1, 0])
```
"""
function getEnergy(mod::CR3BPDynamicsModel, q::Vector{Float64})
    Logging.@debug "getEnergy(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getEnergy" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(q, AbstractVector) && all(isa.(q, Float64)))
        Logging.@error "State vector q must be a Vector{Float64}" type=typeof(q)
        throw(ArgumentError("q must be a Vector{Float64}"))
    end
    if length(q) < 6
        Logging.@error "State vector q is too short" length=length(q) required=6
        throw(ArgumentError("q must have at least 6 elements (position and velocity)"))
    end
    if !all(isfinite, q)
        non_finite_count = count(!isfinite, q)
        Logging.@error "q contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("q must contain only finite values"))
    end

    U::Float64 = getPseudopotential(mod, q[1:3])
    v2::Float64 = q[4]^2 + q[5]^2 + q[6]^2
    if !isfinite(v2) || (v2 < 0)
        Logging.@error "Invalid velocity magnitude" v2=v2
        throw(ArgumentError("Velocity magnitude squared is invalid"))
    end
    JC::Float64 = 2*U-v2
    if !isfinite(JC)
        Logging.@error "Non-finite Jacobi constant computed" JC=JC U=U v2=v2
        throw(ErrorException("Jacobi constant is non-finite"))
    end

    Logging.@info "Computed CR3BP energy" JC=JC pseudopotential=U velocity_sq=v2
    return JC
end

"""
    getEquationsOfMotion(mod::CR3BPDynamicsModel) -> CR3BPEquationsOfMotion

Construct the equations of motion for a CR3BP dynamics model.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `CR3BPEquationsOfMotion`: An equations of motion instance encapsulating the CR3BP differential equations
  for the given dynamics model.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ErrorException` if construction of the equations of motion fails (e.g., invalid primary data).

Notes
- The CR3BP equations describe the motion of a massless test particle in the gravitational field
  of two massive bodies in circular orbits, as viewed from a rotating reference frame.
- All model parameters (masses, scales, etc.) are accessed through the supplied `mod` instance.
- The returned equations of motion object can be used for numerical integration and trajectory analysis.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
eom = getEquationsOfMotion(model)
```
"""
function getEquationsOfMotion(mod::CR3BPDynamicsModel)
    Logging.@debug "getEquationsOfMotion(CR3BP) called"
    
    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getEnergy" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    
    try
        eom = CR3BPEquationsOfMotion(mod)

        Logging.@info "Successfully created CR3BPEquationsOfMotion" primary=mod.primaryData[1].name secondary=mod.primaryData[2].name
        return eom
    catch e
        Logging.@error "Failed to create CR3BPEquationsOfMotion" exception=e
        rethrow()
    end
end

"""
    getEquilibriumPoint(mod::CR3BPDynamicsModel, pointID::Int64) -> Vector{Float64}

Compute the position of a Lagrange equilibrium point in a CR3BP dynamics model rotating frame.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `pointID::Int64`: The equilibrium point identifier (1-5 for L1-L5).

Returns
- `Vector{Float64}`: The 3D position vector of the equilibrium point in normalized coordinates in the rotating frame.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pointID` is not an integer between 1 and 5.
- Throws `ErrorException` if convergence fails for iterative points (L1-L3).
- Throws `ErrorException` if computed position contains non-finite values.

Notes
- Points L1, L2, L3 require iterative solution with Newton-Raphson method.
- Points L4, L5 have analytical solutions.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
L1 = getEquilibriumPoint(model, 1)  # L1 point
```
"""
function getEquilibriumPoint(mod::CR3BPDynamicsModel, pointID::Int64)
    Logging.@debug "getEquilibriumPoint(CR3BP) called" pointID=pointID

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getEquilibriumPoint" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pointID, Integer) && (pointID >= 1) && (pointID <= 5))
        Logging.@error "Invalid equilibrium point ID" pointID=pointID
        throw(ArgumentError("pointID must be an integer between 1 and 5"))
    end

    tol::Float64 = 1E-14
    μ::Float64 = getMassRatios(mod)
    Logging.@debug "Mass ratio for equilibrium computation" μ=μ
    pos::Vector{Float64} = zeros(Float64, 3)
    γ::Float64 = 0.0
    γ_prev::Float64 = -999.0
    count::Int16 = 0
    maxCount::Int16 = 50
    if pointID == 1
        γ = (μ/(3*(1-μ)))^(1/3)
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (μ/γ^2-(1-μ)/(1-γ)^2-γ-μ+1)/(-2*μ/γ^3-2*(1-μ)/(1-γ)^3-1)
            count += 1
        end
        pos[1] = 1-μ-γ
    elseif pointID == 2
        γ = (μ/(3*(1-μ)))^(1/3)
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (-μ/γ^2-(1-μ)/(1+γ)^2-μ+1+γ)/(2*μ/γ^3+2*(1-μ)/(1+γ)^3+1)
            count += 1
        end
        pos[1] = 1-μ+γ
    elseif pointID == 3
        γ = 1-7*μ/12
        while (abs(γ-γ_prev) > tol) && (count < maxCount)
            γ_prev = γ
            γ -= (μ/(-1-γ)^2+(1-μ)/γ^2-μ-γ)/(-2*μ/(1+γ)^3-2*(1-μ)/γ^3-1)
            count += 1
        end
        pos[1] = -μ-γ
    else
        pos[1] = 1/2-μ
        pos[2] = (pointID == 4) ?  sin(pi/3) : -sin(pi/3)
    end
    if (pointID <= 3) && (count >= maxCount)
        Logging.@error "Failed to converge for equilibrium point" pointID=pointID count=count maxCount=maxCount
        throw(ErrorException("Failed to converge to equilibrium point"))
    end
    for i in 1:3
        if !isfinite(pos[i])
            Logging.@error "Non-finite position computed for equilibrium point" pointID=pointID pos=pos
            throw(ErrorException("Equilibrium point position contains non-finite values"))
        end
    end

    Logging.@info "Computed equilibrium point position" pointID=pointID pos=pos
    return pos
end

"""
    getLinearVariationState(mod::CR3BPDynamicsModel, pointID::Int64, var::Vector{Float64}; periodType::String="Short") -> Tuple{Vector{Float64}, Float64}

Compute a linear variation state near an equilibrium point and its associated period in normalized coordinates.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `pointID::Int64`: The equilibrium point identifier (1-5 for L1-L5).
- `var::Vector{Float64}`: A 2- or 3-element position variation vector [Δx, Δy] or [Δx, Δy, Δz].
  If 2 elements provided, z-component is automatically set to 0.
- `periodType::String`: For L4/L5, selects oscillation type: "Short" or "Long" (default: "Short").

Returns
- `Tuple{Vector{Float64}, Float64}`: A tuple containing:
  - The 6D initial state vector at the perturbed location in normalized coordinates.
  - The period of oscillation around the equilibrium point in normalized dimensions.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pointID` is not an integer between 1 and 5.
- Throws `ArgumentError` if `var` is not a `Vector{Float64}` or has fewer than 2 or more than 3 elements.
- Throws `ArgumentError` if variation contains non-finite values.
- Throws `ArgumentError` if `periodType` is not "Short" or "Long".
- Throws `ErrorException` if computed state or period contains non-finite values.

Notes
- For L1-L3 (collinear points): Uses linear approximation for stability analysis.
- For L4-L5 (triangular points): Uses complex eigenvalue method with Short/Long period options.
- Useful for generating initial conditions for periodic orbits and stability studies.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
(q, P) = getLinearVariationState(model, 1, [0.01, 0.01])  # z automatically set to 0
(q, P) = getLinearVariationState(model, 1, [0.01, 0.01, 0])  # explicit z=0
```
"""
function getLinearVariationState(mod::CR3BPDynamicsModel, pointID::Int64, var::Vector{Float64}; periodType::String = "Short")
    Logging.@debug "getLinearVariationState(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getLinearVariationState" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pointID, Integer) && (pointID >= 1) && (pointID <= 5))
        Logging.@error "Invalid equilibrium point ID" pointID=pointID
        throw(ArgumentError("pointID must be an integer between 1 and 5"))
    end
    if !(isa(var, AbstractVector) && all(isa.(var, Float64)))
        Logging.@error "Variation vector var must be a Vector{Float64}" type=typeof(var)
        throw(ArgumentError("var must be a Vector{Float64}"))
    end
    if length(var) < 2 || length(var) > 3
        Logging.@error "Variation vector var has invalid length" length=length(var) valid_range="2-3"
        throw(ArgumentError("var must have 2 or 3 elements"))
    end
    
    # Append 0 for z-component if only 2 elements provided
    if length(var) == 2
        Logging.@debug "Appending z=0 to variation vector" original_length=2
        var = push!(copy(var), 0.0)
    end    
    if !all(isfinite, var)
        non_finite_count = count(!isfinite, var)
        Logging.@error "var contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("var must contain only finite values"))
    end

    μ::Float64 = getMassRatios(mod)
    L::Vector{Float64} = getEquilibriumPoint(mod, pointID)
    ∂²U::Vector{Float64} = getPseudopotentialHessian(mod, L)
    if (1 <= pointID <= 3) # Collinear points L1, L2, L3
        Logging.@debug "Computing linear variation state for collinear point" pointID=pointID
        β_1::Float64 = 2-(∂²U[1]+∂²U[2])/2
        β_2::Float64 = sqrt(-∂²U[1]*∂²U[2])
        s::Float64 = sqrt(β_1+sqrt(β_1^2*β_2^2))
        β_3::Float64 = (s^2+∂²U[1])/(2*s)
        q::Vector{Float64} = [L[1]+var[1], L[2]+var[2], L[3]+var[3], var[2]*s/β_3, -var[1]*s*β_3, 0]
    else  # Triangular points L4, L5
        Logging.@debug "Computing linear variation state for triangular point" pointID=pointID periodType=periodType
        if periodType == "Short"
            λ::Complex{Float64} = sqrt(-1/2-(1/2)*sqrt(complex(1-27*μ*(1-μ))))
            Logging.@debug "Using Short period eigenvalue" λ=λ
        elseif periodType == "Long"
            λ = sqrt(-1/2+(1/2)*sqrt(complex(1-27*μ*(1-μ))))
            Logging.@debug "Using Long period eigenvalue" λ=λ
        else
            Logging.@error "Invalid periodType specified for linear variation state" period=periodType
            throw(ArgumentError("periodType must be 'Short' or 'Long'"))
        end
        s = abs(imag(λ))
        α_1::Float64 = var[1]
        β_1 = var[2]
        α_2::Float64 = (∂²U[4]*α_1+(∂²U[2]+s^2)*β_1)/(2*s)
        β_2 = ((∂²U[1]+s^2)*α_1+∂²U[4]*β_1)/(-2*s)
        q = [L[1]+var[1], L[2]+var[2], L[3]+var[3], α_2*s, β_2*s, 0]
    end
    period::Float64 = 2*π/s
    if !all(isfinite, q)
        Logging.@error "Non-finite state computed" state=q
        throw(ErrorException("Linear variation state contains non-finite values"))
    end
    if !isfinite(period) || period <= 0
        Logging.@error "Invalid period computed" period=period
        throw(ErrorException("Period is non-finite or non-positive"))
    end

    Logging.@info "Computed linear variation state" pointID=pointID periodType=periodType period=period state=q
    return (q, period)
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
    getPrimaryState(mod::CR3BPDynamicsModel, primaryID::Int64) -> Vector{Float64}

Compute the state vector of a primary body in the CR3BP rotating frame.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `primaryID::Int64`: The primary body identifier (1 for primary, 2 for secondary).

Returns
- `Vector{Float64}`: The 6D state vector in normalized rotating coordinates.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `primaryID` is not 1 or 2.
- Throws `ErrorException` if computed state contains non-finite values.

Notes
- In the CR3BP rotating frame, primaries are fixed on the x-axis:
  - Primary 1: at x = -μ
  - Primary 2: at x = (1 - μ)
- Velocity components are all zero since primaries don't move in the rotating frame.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
q_1 = getPrimaryState(model, 1)
```
"""
function getPrimaryState(mod::CR3BPDynamicsModel, primaryID::Int64)
    Logging.@debug "getPrimaryState(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getPrimaryState" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(primaryID, Integer) && (primaryID == 1 || primaryID == 2))
        Logging.@error "Invalid primaryID supplied" primaryID=primaryID
        throw(ArgumentError("primaryID must be 1 or 2"))
    end

    μ::Float64 = getMassRatios(mod)
    q::Vector{Float64} = zeros(Float64, 6)
    q[1] = (primaryID == 1) ? -μ : 1-μ
    if !isfinite(q[1])
        Logging.@error "Non-finite position computed" primaryID=primaryID position=q[1]
        throw(ErrorException("Primary state position is non-finite"))
    end

    Logging.@info "Computed primary state" primaryID=primaryID x_position=q[1] mass_ratio=μ
    return q
end

"""
    getPseudopotential(mod::CR3BPDynamicsModel, pos::Vector{Float64}) -> Float64

Compute the CR3BP pseudopotential (effective potential) at a given position in the rotating frame.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `pos::Vector{Float64}`: A 3D position vector in normalized coordinates.

Returns
- `Float64`: The pseudopotential at the given position.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pos` is not a `Vector{Float64}` or has fewer than 3 elements.
- Throws `ArgumentError` if computed distances are non-finite, zero, or negative.
- Throws `ErrorException` if computed pseudopotential is non-finite.

Notes
- Higher pseudopotential values indicate regions of lower gravitational and
  centrifugal potential.
- Zero crossings of (JC-2*U) define the zero-velocity surfaces (regions of allowed motion).
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
U = getPseudopotential(model, [1.0, 0, 0])
```
"""
function getPseudopotential(mod::CR3BPDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotential(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getPseudopotential" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pos, AbstractVector) && all(isa.(pos, Float64)))
        Logging.@error "Position vector pos must be a Vector{Float64}" type=typeof(pos)
        throw(ArgumentError("pos must be a Vector{Float64}"))
    end
    if length(pos) < 3
        Logging.@error "Position vector pos is too short" length=length(pos) required=3
        throw(ArgumentError("pos must have at least 3 elements (x, y, z coordinates)"))
    end
    if !all(isfinite, pos)
        non_finite_count = count(!isfinite, pos)
        Logging.@error "pos contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("pos must contain only finite values"))
    end

    μ::Float64 = getMassRatios(mod)
    r_13::Float64 = getDistance2Primary(mod, 1, pos)
    r_23::Float64 = getDistance2Primary(mod, 2, pos)
    if !(isfinite(r_13) && (r_13 > 0))
        Logging.@error "Invalid distance to primary" r_13=r_13
        throw(ArgumentError("Distance to primary body is invalid"))
    end
    if !(isfinite(r_23) && (r_23 > 0))
        Logging.@error "Invalid distance to secondary" r_23=r_23
        throw(ArgumentError("Distance to secondary body is invalid"))
    end
    U::Float64 = (1-μ)/r_13+μ/r_23+0.5*(pos[1]^2+pos[2]^2)
    if !isfinite(U)
        Logging.@error "Non-finite pseudopotential computed" U=U r_13=r_13 r_23=r_23 μ=μ pos=pos
        throw(ErrorException("Pseudopotential is non-finite"))
    end

    Logging.@info "Computed CR3BP pseudopotential" U=U r_13=r_13 r_23=r_23
    return U
end

"""
    getPseudopotentialGradient(mod::CR3BPDynamicsModel, pos::Vector{Float64}) -> Vector{Float64}

Compute the gradient (first derivatives) of the CR3BP pseudopotential at a given position.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `pos::Vector{Float64}`: A 3D position vector in normalized rotating coordinates.

Returns
- `Vector{Float64}`: A vector [∂U/∂x, ∂U/∂y, ∂U/∂z] of the pseudopotential
  gradient at the given position.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pos` is not a `Vector{Float64}` or has fewer than 3 elements.
- Throws `ArgumentError` if position contains non-finite values.
- Throws `ErrorException` if computed derivatives contain non-finite values.

Notes
- The gradient represents the force per unit mass in the CR3BP rotating frame.
- Equilibrium points occur where the gradient is zero.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.
- For stability analysis, see `getPseudopotentialJacobian` for second derivatives.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
dU = getPseudopotentialGradient(model, [1.0, 0, 0])
```
"""
function getPseudopotentialGradient(mod::CR3BPDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotentialGradient(CR3BP) called" position=pos

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getPseudopotentialGradient" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pos, AbstractVector) && all(isa.(pos, Float64)))
        Logging.@error "Position vector pos must be a Vector{Float64}" type=typeof(pos)
        throw(ArgumentError("pos must be a Vector{Float64}"))
    end
    if length(pos) < 3
        Logging.@error "Position vector pos is too short" length=length(pos) required=3
        throw(ArgumentError("pos must have at least 3 elements (x, y, z coordinates)"))
    end
    if !all(isfinite, pos)
        non_finite_count = count(!isfinite, pos)
        Logging.@error "pos contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("pos must contain only finite values"))
    end

    μ::Float64 = getMassRatios(mod)
    r_13::Float64 = getDistance2Primary(mod, 1, pos)
    r_23::Float64 = getDistance2Primary(mod, 2, pos)
    if !(isfinite(r_13) && r_13 > 0)
        Logging.@error "Invalid distance to primary 1" r_13=r_13
        throw(ErrorException("Distance to primary 1 is invalid (non-finite or zero)"))
    end
    if !(isfinite(r_23) && r_23 > 0)
        Logging.@error "Invalid distance to primary 2" r_23=r_23
        throw(ErrorException("Distance to primary 2 is invalid (non-finite or zero)"))
    end
    r3_13::Float64 = r_13^3
    r3_23::Float64 = r_23^3
    ∂U∂x::Float64 = pos[1]-(1-μ)*(pos[1]+μ)/r3_13-μ*(pos[1]+μ-1)/r3_23
    ∂U∂y::Float64 = pos[2]-(1-μ)*pos[2]/r3_13-μ*pos[2]/r3_23
    ∂U∂z::Float64 = -(1-μ)*pos[3]/r3_13-μ*pos[3]/r3_23
    if !isfinite(∂U∂x)
        Logging.@error "Non-finite derivative computed" component="x" value=∂U∂x
        throw(ErrorException("∂U/∂x is non-finite"))
    end
    if !isfinite(∂U∂y)
        Logging.@error "Non-finite derivative computed" component="y" value=∂U∂y
        throw(ErrorException("∂U/∂y is non-finite"))
    end
    if !isfinite(∂U∂z)
        Logging.@error "Non-finite derivative computed" component="z" value=∂U∂z
        throw(ErrorException("∂U/∂z is non-finite"))
    end

    Logging.@info "Computed CR3BP pseudopotential gradient" ∂U∂x=∂U∂x ∂U∂y=∂U∂y ∂U∂z=∂U∂z
    return [∂U∂x, ∂U∂y, ∂U∂z]
end

"""
    getPseudopotentialHessian(mod::CR3BPDynamicsModel, pos::Vector{Float64}) -> Vector{Float64}

Compute the Hessian (second derivatives) of the CR3BP pseudopotential at a given position.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.
- `pos::Vector{Float64}`: A 3D position vector in normalized rotating coordinates.

Returns
- `Vector{Float64}`: A vector [∂²U/∂x², ∂²U/∂y², ∂²U/∂z², ∂²U/∂x∂y, ∂²U/∂x∂z, ∂²U/∂y∂z]
  containing the Hessian matrix elements of the pseudopotential in order.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if `pos` is not a `Vector{Float64}` or has fewer than 3 elements.
- Throws `ArgumentError` if position contains non-finite values.
- Throws `ErrorException` if computed second derivatives contain non-finite values.

Notes
- The Hessian is essential for stability analysis of equilibrium points and orbits.
- Eigenvalues of the Hessian matrix determine linear stability characteristics.
- At equilibrium points (where first derivatives are zero), the Hessian defines
  the center manifold structure.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.
- Use `getPseudopotentialDerivatives` to compute first derivatives.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
ddU = getPseudopotentialHessian(model, [1.0, 0, 0])
```
"""
function getPseudopotentialHessian(mod::CR3BPDynamicsModel, pos::Vector{Float64})
    Logging.@debug "getPseudopotentialHessian(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to getPseudopotentialHessian" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end
    if !(isa(pos, AbstractVector) && all(isa.(pos, Float64)))
        Logging.@error "Position vector pos must be a Vector{Float64}" type=typeof(pos)
        throw(ArgumentError("pos must be a Vector{Float64}"))
    end
    if length(pos) < 3
        Logging.@error "Position vector pos is too short" length=length(pos) required=3
        throw(ArgumentError("pos must have at least 3 elements (x, y, z coordinates)"))
    end
    if !all(isfinite, pos)
        non_finite_count = count(!isfinite, pos)
        Logging.@error "pos contains non-finite values" non_finite_count=non_finite_count
        throw(ArgumentError("pos must contain only finite values"))
    end

    μ::Float64 = getMassRatios(mod)
    r_13::Float64 = getDistance2Primary(mod, 1, pos)
    r_23::Float64 = getDistance2Primary(mod, 2, pos)
    if !(isfinite(r_13) && r_13 > 0)
        Logging.@error "Invalid distance to primary 1" r_13=r_13
        throw(ErrorException("Distance to primary 1 is invalid (non-finite or zero)"))
    end
    if !(isfinite(r_23) && r_23 > 0)
        Logging.@error "Invalid distance to primary 2" r_23=r_23
        throw(ErrorException("Distance to primary 2 is invalid (non-finite or zero)"))
    end
    r3_13::Float64 = r_13^3
    r3_23::Float64 = r_23^3
    r5_13::Float64 = r_13^5
    r5_23::Float64 = r_23^5
    ∂²U∂x²::Float64 = 1-(1-μ)/r3_13-μ/r3_23+3*(1-μ)*(pos[1]+μ)^2/r5_13+3*μ*(pos[1]+μ-1)^2/r5_23
    ∂²U∂y²::Float64 = 1-(1-μ)/r3_13-μ/r3_23+3*(1-μ)*pos[2]^2/r5_13+3*μ*pos[2]^2/r5_23
    ∂²U∂z²::Float64 = -(1-μ)/r3_13-μ/r3_23+3*(1-μ)*pos[3]^2/r5_13+3*μ*pos[3]^2/r5_23
    ∂²U∂x∂y::Float64 = 3*(1-μ)*(pos[1]+μ)*pos[2]/r5_13+3*μ*(pos[1]+μ-1)*pos[2]/r5_23
    ∂²U∂x∂z::Float64 = 3*(1-μ)*(pos[1]+μ)*pos[3]/r5_13+3*μ*(pos[1]+μ-1)*pos[3]/r5_23
    ∂²U∂y∂z::Float64 = 3*(1-μ)*pos[2]*pos[3]/r5_13+3*μ*pos[2]*pos[3]/r5_23
    ∂²U::Vector{Float64} = [∂²U∂x², ∂²U∂y², ∂²U∂z², ∂²U∂x∂y, ∂²U∂x∂z, ∂²U∂y∂z]
    ∂²U_names::Vector{String} = ["∂²U/∂x²", "∂²U/∂y²", "∂²U/∂z²", "∂²U/∂x∂y", "∂²U/∂x∂z", "∂²U/∂y∂z"]
    for (i, (val, name)) in enumerate(zip(∂²U, ∂²U_names))
        if !isfinite(val)
            Logging.@error "Non-finite second derivative computed" component=name value=val
            throw(ErrorException("$name is non-finite"))
        end
    end

    Logging.@info "Computed CR3BP pseudopotential Jacobian" ∂²U∂x²=∂²U[1] ∂²U∂y²=∂²U[2] ∂²U∂z²=∂²U[3]
    return ∂²U
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
    - `FULL → 42`
    - `ARCLENGTH → 43`
    - `MOMENTUM → 43`
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

    sizeMap = Dict(SIMPLE => 6, STM => 42, FULL => 42, ARCLENGTH => 43, MOMENTUM => 43)
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
    isEpochIndependent(mod::CR3BPDynamicsModel) -> Bool

Determine that a CR3BP dynamics model's equations of motion are epoch-independent.

Arguments
- `mod::CR3BPDynamicsModel`: A valid CR3BP dynamics model with exactly 2 primaries.

Returns
- `Bool`: Always returns `true` for CR3BP, since the rotating frame formulation is time-autonomous.

Errors
- Throws `ArgumentError` if `mod` is not a valid `CR3BPDynamicsModel` or does not
  contain exactly 2 primary bodies.
- Throws `ArgumentError` if the model structure is invalid.

Notes
- The CR3BP in the rotating frame is time-invariant (autonomous).
- This means the equations of motion do not explicitly depend on the epoch/time.
- Propagation can use fixed time-step methods independent of the absolute epoch.
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
model = CR3BPDynamicsModel(sys, [1, 2])
indep = isEpochIndependent(model)  # Returns true
```
"""
function isEpochIndependent(mod::CR3BPDynamicsModel)::Bool
    Logging.@debug "isEpochIndependent(CR3BP) called"

    if mod === nothing || !isa(mod, CR3BPDynamicsModel)
        Logging.@error "Invalid model supplied to isEpochIndependent" type=typeof(mod)
        throw(ArgumentError("mod must be a CR3BPDynamicsModel instance"))
    end
    if length(mod.primaryData) != 2
        Logging.@error "CR3BP requires exactly two primaries" n=length(mod.primaryData)
        throw(ArgumentError("CR3BPDynamicsModel must contain exactly two primary bodies"))
    end

    # CR3BP is epoch-independent (autonomous in rotating frame)
    is_independent::Bool = true

    Logging.@info "CR3BPDynamicsModel is epoch independent" result=is_independent
    return is_independent
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
