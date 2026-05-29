"""
Dynamics model types

Author: Jonathan LeFevre Richmond
C: 5/1/26
U: 5/29/26
"""


"""
Abstract type for dynamics model types
"""
abstract type AbstractDynamicsModel end


"""
    init_dynamicsModel(systemData::SystemData, indices::Vector{Int64}, ::Type{M}) -> M <: AbstractDynamicsModel

Initialize dynamics model from given list of bodies in a system data object

# Arguments
- `systemData::SystemData`: System data object
- `indices::Vector{Int64}`: List of indices specifying which bodies to include
- `::Type{M}`: Concrete dynamics model type

Returns
- `M<:AbstractDynamicsModel`: Concrete dynamics model object

Errors
- Throws `ArgumentError` if `indices` is empty or if there are duplicate
    indices
- Throws `BoundError` if element of `indices` is out of range of the length of
    `systemData.bodyData`

Logging
- Emits `@error` logs for thrown errors
- Emits `@info` logs when dynamics model is initialized
- Emits `@debug` logs when function is entered, when initializing dynamics
    model, or when primary bodies are extracted

Example
```
dynamicsModel = init_dynamicsModel(systemData, [1, 2], CR3BPDynamicsModel)
```
"""
function init_dynamicsModel(systemData::SystemData, indices::Vector{Int64}, ::Type{M})::M where {M <: AbstractDynamicsModel}
    Logging.@debug "Entered init_dynamicsModel" systemData indices modelType=M
    
    # Input validation
    if isempty(indices)
        Logging.@error "Index list must be non-empty" modelType=M
        throw(ArgumentError("Indices cannot be empty"))
    end

    n::Int64 = length(systemData.bodyData)
    for idx in indices
        if (idx < 1) || (idx > n)
            Logging.@error "Body index out of range" idx validRange=(1, n) modelType=M
            throw(BoundsError(systemData.bodyData, idx))
        end
    end
    if length(unique(indices)) != length(indices)
        Logging.@error "Duplicate indices detected" indices
        throw(ArgumentError("Indices contains duplicates"))
    end

    Logging.@debug "Initializing dynamics model" modelType=M indices totalBodies=n

    # Extract body data (pre-allocated, single pass)
    primaryData = Vector{BodyData}(undef, length(indices))
    for (i::Int64, idx::Int64) in enumerate(indices)
        primaryData[i] = systemData.bodyData[idx]
    end

    Logging.@debug "Extracted primary bodies" names=[b.name for b in primaryData] modelType=M

    # Dispatch to model-specific constructor
    local model::M
    try
        model = build_dynamicsModel(M, primaryData)
    catch e
        Logging.@error "Failed to construct dynamics model" modelType=M exception=(e, catch_backtrace())
        rethrow()
    end

    Logging.@info "Dynamics model initialized" modelType=M bodies=[b.name for b in primaryData]

    return model
end


"""
    CR3BPDynamicsModel(systemData::SystemData, indices::Vector{Int64})

Fields
- `primaryData::Vector{BodyData}`: List of primary body data

Construction
- Direct: `CR3BPDynamicsModel(primaryData)`
- Convenience: `CR3BPDynamicsModel(systemData, indices)` initializes CR3BP
    dynamics model via `init_dynamicsModel`

Example
```
dynamicsModel = CR3BPDynamicsModel(systemData, [1, 2])
println(dynamicsModel)
```
"""
struct CR3BPDynamicsModel <: AbstractDynamicsModel
    primaryData::Vector{BodyData}

    function CR3BPDynamicsModel(primaryData::Vector{BodyData})
        if length(primaryData) != 2
            Logging.@error "CR3BP requires exactly 2 primary bodies" n=length(primaryData)
            throw(ArgumentError("CR3BPDynamicsModel requires exactly 2 bodies, got $(length(primaryData))"))
        end
        primary::BodyData = primaryData[1]
        secondary::BodyData = primaryData[2]
        if secondary.parentSpiceID != primary.spiceID
            Logging.@error "CR3BP requires primary to be parent of secondary" primary=primary.name primarySpice=primary.spiceID secondary=secondary.name secondaryParent=secondary.parentSpiceID
            throw(ArgumentError("CR3BPDynamicsModel requires secondary's parent SPICE ID to equal primary's SPICE ID"))
        end

        new(primaryData)
    end
end
CR3BPDynamicsModel(systemData::SystemData, indices::Vector{Int64}) = init_dynamicsModel(systemData, indices, CR3BPDynamicsModel)

"""
    build_dynamicsModel(::Type{CR3BPDynamicsModel}, primaryData::Vector{BodyData}) -> CR3BPDynamicsModel

Build CR3BP dynamics model from given list of bodies

Arguments
- `::Type{CR3BPDynamicsModel}`: CR3BPDynamicsModel type
- `primaryData::Vector{BodyData}`: Bodies

Returns
- `CR3BPDynamicsModel`: CR3BPDynamicsModel object

Errors
- Throws `ArgumentError` if there are not exactly 2 bodies in `primaryData` or
    if the primary body is not the parent of the secondary

Logging
- Emits `@error` logs for thrown errors
- Emits `@debug` logs when function is entered

Example
```
dynamicsModel_CR3BP = build_dynamicsModel(CR3BPDynamicsModel, [Earth, Moon])
```
"""
function build_dynamicsModel(::Type{CR3BPDynamicsModel}, primaryData::Vector{BodyData})::CR3BPDynamicsModel
    Logging.@debug "Entered build_dynamicsModel for CR3BPDynamicsModel" primaryData
    
    # CR3BP requires exactly two bodies
    if length(primaryData) != 2
        Logging.@error "CR3BP requires exactly 2 primary bodies" n=length(primaryData)
        throw(ArgumentError("CR3BPDynamicsModel requires exactly 2 bodies, got $(length(primaryData))"))
    end

    # Secondary's parent must be the primary
    primary::BodyData = primaryData[1]
    secondary::BodyData = primaryData[2]
    if secondary.parentSpiceID != primary.spiceID
        Logging.@error "CR3BP requires primary to be parent of secondary" primary=primary.name primarySpice=primary.spiceID secondary=secondary.name secondaryParent=secondary.parentSpiceID
        throw(ArgumentError("CR3BPDynamicsModel requires secondary's parent SPICE ID to equal primary's SPICE ID"))
    end

    return CR3BPDynamicsModel(primaryData)
end


# Fixed mapping from EquationType to state vector size
const _CR3BP_state_sizes = Dict{EquationType, Int64}(
    SIMPLE      => 6,
    STM         => 42,
    ARCLENGTH   => 7,
    MOMENTUM    => 7,
    FULL        => 44
)


# Base functions
Base.:(==)(dynamicsModel1::AbstractDynamicsModel, dynamicsModel2::AbstractDynamicsModel) = (dynamicsModel1.primaryData == dynamicsModel2.primaryData)
Base.isequal(dynamicsModel1::AbstractDynamicsModel, dynamicsModel2::AbstractDynamicsModel) = isequal(dynamicsModel1.primaryData, dynamicsModel2.primaryData)
function Base.show(io::IO, ::MIME"text/plain", dynamicsModel::CR3BPDynamicsModel)
    n = length(dynamicsModel.primaryData)
    println(io, "CR3BPDynamicsModel: ", n, n == 1 ? " primary" : " primaries")
    for body in dynamicsModel.primaryData
        println(io, "  ─────────────────────────────────────")
        println(io, "  Body: ", body.name)
        println(io, "    SPICE ID: ", body.spiceID, "   Parent SPICE ID: ", body.parentSpiceID)
        Printf.@printf(io, "    a: %0.6g km   e: %0.6g   i: %0.6g rad\n", body.a, body.e, body.i)
        Printf.@printf(io, "    radius: %0.6g km   μ (GM): %0.6g km^3/s^2   mass: %0.6g kg\n", body.r, body.μ, body.m)
        Printf.@printf(io, "    Ω (RAAN): %0.6g rad\n", body.Ω)
    end
    println(io, "  ─────────────────────────────────────")
end
function Base.show(io::IO, dynamicsModel::AbstractDynamicsModel)
    Base.show(io, MIME"text/plain"(), dynamicsModel)
end
