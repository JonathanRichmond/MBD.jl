"""
SPICE utility functions

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 4/16/26
"""


# Cache for memoization
const _id_cache = Dict{String, Int64}()


"""
    getIDCode(name::String) -> Int64

Get the SPICE ID for a given body name

Arguments
- `name::String`: Body name (as accepted by SPICE)

Returns
- `Int64`: SPICE ID corresponding to `name` as returned by `SPICE.bods2c`

Errors
- Throws `ArgumentError` if `name` is empty
- Throws `ErrorException` when `SPICE.bods2c` call fails
- Throws `KeyError` if no SPICE ID is returned

Logging
- Emits `@error` logs for thrown errors
- Emits `@warn` logs if no SPICE ID is returned
- Emits `@info` logs when a SPICE ID is resolved
- Emits `@debug` logs when function is entered, if `name` is cached, or when
    SPICE is queried

Notes
- Function has an injectable reference for testing

Example
```
id = getIDCode("Earth")
```
"""
function getIDCode_func(name::String)
    Logging.@debug "Entered getIDCode_func()" name

    # Input validation
    if isempty(strip(name))
        Logging.@error "Body name must be a non-empty string" name
        throw(ArgumentError("Body name cannot be empty or whitespace"))
    end

    # Normalize name once to avoid redundant allocations
    normalized = String(strip(name))

    # Cache lookup (avoids redundant SPICE calls)
    if haskey(_id_cache, normalized)
        Logging.@debug "Cache hit for body name" name=normalized id=_id_cache[normalized]
        return _id_cache[normalized]
    end

    Logging.@debug "Cache miss - querying SPICE" name=normalized

    # SPICE lookup
    local code
    try
        code = SPICE.bods2c(normalized)
    catch e
        Logging.@error "SPICE lookup failed for body name" name=normalized exception=(e, catch_backtrace())
        rethrow()
    end

    # Error if name is not found
    if isnothing(code)
        Logging.@warn "No SPICE Id found for body name" name=normalized
        throw(KeyError("No SPICE ID found for body: \"$normalized\""))
    end

    # Store result in cache
    _id_cache[normalized] = code
    Logging.@info "Resolved SPICE ID" name=normalized id=code

    return code
end
# Injectable function reference
const _getIDCode_func = Ref{Function}(getIDCode_func)
# Callable function
function getIDCode(name::String)
    code = _getIDCode_func[](name)

    return code
end
