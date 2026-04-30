"""
SystemData methods

Author: Jonathan LeFevre Richmond
C: 4/30/26
"""


"""
    getNumPrimaries(systemData::SystemData) -> Int64

Return number of primary bodies contained in a `SystemData` object

Arguments
- `systemData::SystemData`: `SystemData` object

Returns
- `Int64`: Number of entries in `systemData.bodyData`

Errors
- Throws `ArgumentError` if `systemData` has inconsistent field lengths

Logging
- Emits `@error` logs for thrown errors
- Emits `@warn` logs if `systemData.bodyData` is empty
- Emits `@debug` logs when function is entered or when `systemData.bodyData`
    contains entries

Example
```
nPrimaries::Int64 = getNumPrimaries(systemData)
```
"""
function getNumPrimaries(systemData::SystemData)::Int64
    Logging.@debug "Entered getNumPrimaries" systemData

    n::Int64 = length(systemData.bodyData)
    if n == 0
        Logging.@warn "SystemData contains no primaries (bodyData is empty)"
    else
        Logging.@debug "SystemData contains primaries" count=n
    end

    if (length(systemData.names) != n) || (length(systemData.spiceIDs) != n)
        Logging.@error "SystemData has inconsistent field lengths" bodyData_len=n names_len=length(systemData.names) spiceIDs_len=length(systemData.spiceIDs)
        throw(ArgumentError("SystemData fields bodyData, names, and spiceIDs must have equal length"))
    end
    
    return n
end


"""
    shallowClone(systemData::SystemData) -> SystemData

Return shallow copy of `SystemData` object

Arguments
- `systemData::SystemData`: `SystemData` object

Returns
- `SystemData`: Copy of `SystemData` object

Errors
- No additional error checking

Logging
- Emits `@debug` logs when function is entered or when shallow copy is created

Example
```
systemData2::SystemData = shallowClone(systemData)
```
"""
function shallowClone(systemData::SystemData)::SystemData
    Logging.@debug "Entered shallowClone" systemData

    clone = SystemData(copy(systemData.bodyData), copy(systemData.names), copy(systemData.spiceIDs))

    Logging.@debug "Shallow copy created" nPrimaries=length(clone.bodyData)

    return clone
end
