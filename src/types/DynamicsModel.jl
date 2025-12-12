"""
Dynamics model types

Author: Jonathan Richmond
C: 12/12/25
"""


"""
    init_dynamicsModel(systemData::SystemData, indices::Vector{Int64}, modelType::ModelType)

Initialize a dynamics model for a subset of bodies from a system.

Arguments
- `systemData::SystemData`: The system data containing all available bodies.
- `indices::Vector{Int64}`: 1-based indices into `systemData.bodyData` specifying
  which bodies to include as primaries in the dynamics model. Must be non-empty,
  within bounds, and contain no duplicates.
- `modelType::ModelType`: The type of dynamics model to construct (e.g., `CR3BP`).

Returns
- A dynamics model instance appropriate for the specified `modelType`. For `CR3BP`,
  returns a `CR3BPDynamicsModel` constructed from the primary bodies.

Errors
- Throws `ArgumentError` if:
  - `systemData` is not a valid `SystemData` instance
  - `indices` is not a `Vector{Int64}` or is empty
  - `indices` contains duplicates
  - Model-specific requirements are not met (e.g., CR3BP requires exactly 2 bodies
    where the secondary's parent SPICE ID matches the primary's SPICE ID)
  - `modelType` is unsupported
- Throws `BoundsError` if any index in `indices` is outside [1, length(systemData.bodyData)]

Model-specific requirements
- `CR3BP`: Requires exactly 2 indices. The first is the primary body (e.g., planet),
  the second is the secondary (e.g., moon). The secondary's `parentSPICEID` must
  equal the primary's `SPICEID`.

Notes
- The function emits `Logging` messages at `@debug`, `@info`, and `@error` levels
  to aid troubleshooting.

Example
```
sys = SystemData(["Earth", "Moon"])
model = init_dynamicsModel(sys, [1, 2], CR3BP)
```
"""
function init_dynamicsModel(systemData::MBD.SystemData, indices::Vector{Int64}, modelType::MBD.ModelType)
    Logging.@debug "init_dynamicsModel called" n_indices=length(indices) modelType=modelType

    # Basic validation
    if systemData === nothing || !isa(systemData, MBD.SystemData)
        Logging.@error "Invalid systemData supplied to init_dynamicsModel"
        throw(ArgumentError("systemData must be a MBD.SystemData instance"))
    end
    if !(isa(indices, AbstractVector) && all(isa.(indices, Int64)))
        Logging.@error "indices must be a Vector{Int64}" indices=indices
        throw(ArgumentError("indices must be a Vector{Int64}"))
    end
    if isempty(indices)
        Logging.@error "indices is empty" 
        throw(ArgumentError("indices must be non-empty"))
    end

    nBodies = length(systemData.bodyData)
    # Validate index bounds and uniqueness
    if any(x -> (x < 1) || (x > nBodies), indices)
        Logging.@error "One or more indices out of bounds" indices=indices nBodies=nBodies
        throw(BoundsError("indices must be between 1 and $(nBodies)"))
    end
    if length(unique(indices)) != length(indices)
        Logging.@error "Duplicate indices supplied" indices=indices
        throw(ArgumentError("indices must not contain duplicates"))
    end

    # Collect primary bodies
    primaryData::Vector{MBD.BodyData} = []
    for idx in indices
        push!(primaryData, systemData.bodyData[idx])
    end

    # Dispatch by model type with model-specific checks
    if modelType == CR3BP
        # CR3BP requires exactly two bodies: primary then secondary
        if length(primaryData) != 2
            Logging.@error "CR3BP requires exactly two primary bodies" provided=length(primaryData)
            throw(ArgumentError("CR3BP model requires exactly two bodies (primary, secondary)"))
        end
        # Secondary's parent must be the primary
        primary = primaryData[1]
        secondary = primaryData[2]
        if secondary.parentSPICEID != primary.SPICEID
            Logging.@error "CR3BP primary/secondary parent mismatch" primary=primary.name primarySPICE=primary.SPICEID secondary=secondary.name secondaryParent=secondary.parentSPICEID
            throw(ArgumentError("CR3BP requires the secondary's parent SPICE ID to equal the primary's SPICE ID"))
        end
        Logging.@info "Initializing CR3BP dynamics model" primary=primary.name secondary=secondary.name
        return CR3BPDynamicsModel(primaryData)
    end

    # Unknown model type
    Logging.@error "Unsupported model type" modelType=modelType
    throw(ArgumentError("unsupported model type: $(modelType)"))
end


abstract type AbstractDynamicsModel end


"""
    CR3BPDynamicsModel <: AbstractDynamicsModel

Dynamics model for the Circular Restricted 3-Body Problem (CR3BP).

Fields
- `primaryData::Vector{BodyData}`: Two-element vector where index 1 is the primary
  body (e.g., planet) and index 2 is the secondary body (e.g., moon). The
  secondary's `parentSPICEID` must equal the primary's `SPICEID`.

Construction
- Direct construction: `CR3BPDynamicsModel(primaryData::Vector{BodyData})`
- Convenience constructor: `CR3BPDynamicsModel(systemData::SystemData, indices::Vector{Int64})`
  calls `init_dynamicsModel(systemData, indices, CR3BP)` to validate inputs and
  construct the model. The `indices` argument must contain exactly two unique,
  in-bounds indices referencing bodies in `systemData.bodyData`.

Notes
- The CR3BP model assumes circular orbits and restricts motion to the orbital plane
  of the two primary bodies.
- Proper parent-child relationship validation ensures the secondary orbits the primary.

Example
```
sys = SystemData(["Earth", "Moon"])
model = CR3BPDynamicsModel(sys, [1, 2])
println(model)  # displays formatted model information
```
"""
struct CR3BPDynamicsModel <: AbstractDynamicsModel
    primaryData::Vector{MBD.BodyData}
end
CR3BPDynamicsModel(systemData::MBD.SystemData, indices::Vector{Int64}) = init_dynamicsModel(systemData, indices, MBD.CR3BP)
Base.:(==)(mod1::MBD.CR3BPDynamicsModel, mod2::MBD.CR3BPDynamicsModel) = (mod1.primaryData == mod2.primaryData)
function Base.show(io::IO, ::MIME"text/plain", mod::MBD.CR3BPDynamicsModel)
    println(io, "CR3BPDynamicsModel:")
    n = length(mod.primaryData)
    println(io, "  Primaries: ", n)
    if n >= 1
        primary = mod.primaryData[1]
        Printf.@printf(io, "    1) %-19s  SPICEID=%6d\n", primary.name, primary.SPICEID)
    end
    if n >= 2
        secondary = mod.primaryData[2]
        Printf.@printf(io, "    2) %-19s  SPICEID=%6d  parent=%6d\n", secondary.name, secondary.SPICEID, secondary.parentSPICEID)
    end
end
function Base.show(io::IO, mod::MBD.CR3BPDynamicsModel)
    Base.show(io, MIME"text/plain"(), mod)
end
