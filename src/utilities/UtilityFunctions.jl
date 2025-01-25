"""
Utility functions

Author: Jonathan Richmond
C: 9/7/22
U: 1/25/25
"""

export Cartesian2Cylindrical, checkIndices, isApproxSigFigs, maskData, updatePointer

"""
    Cartesian2Cylindrical(pos, c)

Return position relative to point in cylindrical coordinates

# Arguments
- `pos::Vector{Float64}`: Position vector
- `c::Vector{Float64}`: Center of cylinder
"""
function Cartesian2Cylindrical(pos::Vector{Float64}, c::Vector{Float64})
    return [sqrt((pos[1]-c[1])^2+(pos[2]-c[2])^2), atan(pos[2]-c[2], pos[1]-c[1]), pos[3]-c[3]]
end

"""
    checkIndices(stateIndices, stateSize)

Throw error if indices are invalid

# Arguments
- `stateIndices::Vector{Int64}`: Indices of state vector
- `stateSize::Int64`: Number of state get_elements_by_tagname
"""
function checkIndices(stateIndices::Vector{Int64}, stateSize::Int64)
    state::StaticArrays.SVector{stateSize, Float64} = StaticArrays.SVector{stateSize, Float64}(zeros(Float64, stateSize))
    for i::Int16 in Int16(1):Int16(length(stateIndices))
        ((stateIndices[i] < 1) || (stateIndices[i] > stateSize)) && throw(BoundsError(state, stateIndices[i]))
        for j::Int16 in Int16(1):Int16(length(stateIndices))
            ((i != j) && (stateIndices[i] == stateIndices[j])) && throw(ArgumentError("Index $(stateIndices[i]) is included more than once"))
        end
    end
end

"""
    isApproxSigFigs(a, b, sigFigs)

# Arguments
- `a::Vector{Float64}`: First number
- `b::Vector{Float64}`: Second number
- `sigFigs::Int64`: Number of significant figures
"""
function isApproxSigFigs(a::Vector{Float64}, b::Vector{Float64}, sigFigs::Int64)
    (length(a) == length(b)) || throw(ArgumentError("Length of a, $(length(a)), must match lenght of b, $(length(b))"))
    for i::Int16 = Int16(1):Int16(length(a))
        (round(a[i], RoundNearestTiesUp, sigdigits = sigFigs) == round(b[i], RoundNearestTiesUp, sigdigits = sigFigs)) || (return false)
    end

    return true
end

"""
    maskData(mask, data)

Return masked data

# Arguments
- `mask::Vector{Bool}`: Column mask
- `data::Matrix{Float64}`: Data
"""
function maskData(mask::Vector{Bool}, data::Matrix{Float64})
    (isempty(mask) || (size(data, 2) == 0)) && (return zeros(Float64, (1, 0)))
    n_col::Int16 = Int16(length(filter(x -> x == true, mask)))
    (size(data, 2) == length(mask)) || throw(ArgumentError("Mask length, $(length(mask)), must match data row length, $(size(data, 2))"))
    if n_col == Int16(length(mask))
        return data
    elseif n_col == Int16(0)
        return zeros(Float64, (1, 0))
    else
        count::Int16 = 1
        maskedData::Matrix{Float64} = zeros(Float64, (size(data, 1),n_col))
        for c::Int16 in Int16(1):Int16(length(mask))
            if mask[c]
                [(maskedData[r,count] = data[r,c]) for r in 1:size(data, 1)]
                count += 1
            end
        end

        return maskedData
    end
end

"""
    updatePointer(original, copiedObjectMap, forceMatch)

Update pointer

# Arguments
- `original::Any`: Object
- `copiedObjectMap::Dict`: Map between old and new objects
- `forceMatch::Bool`: Force match?
"""
function updatePointer(original::Any, copiedObjectMap::Dict, forceMatch::Bool)
    contains::Bool = false
    for key in keys(copiedObjectMap)
        if key == hash(original)
            contains = true
            break
        end
    end
    if contains
        return copiedObjectMap[hash(original)]
    else
        forceMatch ? throw(ErrorException("Could not find match for original in copiedObjectMap")) : (return original)
    end
end
