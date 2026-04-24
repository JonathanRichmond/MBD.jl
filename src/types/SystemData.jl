"""
System data types

Author: Jonathan LeFevre Richmond
C: 4/24/26
"""


"""
    init_systemData(names::Vector{String}) -> SystemData

Initialize system data from a given list of body names

Arguments
- `names::Vector{String}`: Body names (as accepted by SPICE)

Returns
- `SystemData`: SystemData struct

Errors
- Throws `ArgumentError` if `names` is empty, if an element of `names` is
    empty, or if there are duplicate SPICE IDs

Logging
- Emits `@error` logs for thrown errors
- Emits `@warn` logs if duplicate body names are detected
- Emits `@info` logs when system data is initialized
- Emits `@debug` logs when function is entered, when iniitalizing system data,
    when loading body data, or when body data is loaded

Example
```
systemData = init_systemData(["Earth", "Moon"])
```
"""
function init_systemData(names::Vector{String})::SystemData
    Logging.@debug "Entered init_systemData" names

    # Input validation
    if isempty(names)
        Logging.@error "Body name list must be non-empty"
        throw(ArgumentError("Body names vector cannot be empty"))
    end

    # Validate and normalize all names before any loading
    normalized = Vector{String}(undef, length(names))
    for (idx::Int64, name::String) in enumerate(names)
        if isempty(strip(name))
            Logging.@error "Body name at index $idx must be a non-empty string" name
            throw(ArgumentError("Body name at index $idx cannot be empty or whitespace"))
        end
        normalized[idx] = String(strip(name))
    end

    # Warn about duplicates
    uniqueNames::Vector{String} = unique(normalized)
    if length(uniqueNames) < length(normalized)
        duplicates::Vector{String} = filter(n -> count(==(n), normalized) > 1, uniqueNames)
        Logging.@warn "Duplicate body names detected" duplicates
    end

    Logging.@debug "Initializing system data" bodies=normalized total=length(normalized)

    # Preallocate with known size to avoid repeated heap allocations
    nNames::Int64 = length(normalized)
    bodyData = Vector{BodyData}(undef, nNames)
    spiceIDs = Vector{Int64}(undef, nNames)

    # Load body data for each name
    for (idx::Int64, name::String) in enumerate(normalized)
        Logging.@debug "Loading body data" name index=idx total=nNames
        local data::BodyData
        try
            data = BodyData(name)
        catch e
            Logging.@error "Failed to load body data" name index=idx exception=(e, catch_backtrace())
            rethrow()
        end
        bodyData[idx] = data
        spiceIDs[idx] = data.spiceID
        Logging.@debug "Loaded body data" name SPICEID=data.spiceID index=idx
    end

    # Validate that all SPICE IDs are unique
    if length(unique(spiceIDs)) < nNames
        Logging.@error "Duplicate SPICE IDs detected across bodies" spiceIDs names=normalized
        throw(ArgumentError("Loaded bodies contain duplicate SPICE IDs: $spiceIDs"))
    end

    result = SystemData(bodyData, normalized, spiceIDs)
    Logging.@info "System data initialized" bodies=normalized spiceIDs total=nNames

    return result
end


"""
    SystemData(names::Vector{String})

Fields
- `bodyData::Vector{BodyData}`: List of body data corresponding to `names`
- `names::Vector{String}`: Body names (as accepted by SPICE)
- `spiceIDs::Vector{Int64}`: List of SPICE IDs corresponding to `names`

Construction
- Direct: `SystemData(bodyData, names, spiceIDs)`
- Convenience: `SystemData(names)` initializes system data via
    `init_systemData`

Example
```
systemData = SystemData(["Earth", "Moon"])
println(systemData)
```
"""
struct SystemData
    bodyData::Vector{BodyData}
    names::Vector{String}
    spiceIDs::Vector{Int64}
end
SystemData(names::Vector{String}) = init_systemData(names)


# Base functions
Base.:(==)(systemData1::SystemData, systemData2::SystemData) = (systemData1.bodyData == systemData2.bodyData) && (systemData1.names == systemData2.names) && (systemData1.spiceIDs == systemData2.spiceIDs)
function Base.show(io::IO, ::MIME"text/plain", systemData::SystemData)
    nNames::Int64 = length(systemData.names)
    println(io, "SystemData: ", nNames, nNames == 1 ? " body" : " bodies")
    for (name, id, body) in zip(systemData.names, systemData.spiceIDs, systemData.bodyData)
        println(io, "  ─────────────────────────────────────")
        println(io, "  Body: ", name)
        println(io, "    SPICE ID: ", id, "   Parent SPICE ID: ", body.parentSpiceID)
        Printf.@printf(io, "    a: %0.6g km   e: %0.6g   i: %0.6g rad\n", body.a, body.e, body.i)
        Printf.@printf(io, "    radius: %0.6g km   μ (GM): %0.6g km^3/s^2   mass: %0.6g kg\n", body.r, body.μ, body.m)
        Printf.@printf(io, "    Ω (RAAN): %0.6g rad\n", body.Ω)
    end
    println(io, "  ─────────────────────────────────────")
end
function Base.show(io::IO, systemData::SystemData)
    Base.show(io, MIME"text/plain"(), systemData)
end
