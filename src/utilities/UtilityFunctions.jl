"""
Utility functions

Author: Jonathan Richmond
C: 9/7/22
U: 4/3/24
"""

export Cartesian2Cylindrical, checkIndices, maskData, updatePointer

"""
    Cartesian2Cylindrical(pos, c)

Return position relative to point in cylindrical coordinates

# Arguments
- `pos::Vector{Float64}`: Position vector
- `c::Vector{Float64}`: Center of cylinder
"""
function Cartesian2Cylindrical(pos::Vector{Float64}, c::Vector{Float64})
    r::Float64 = sqrt((pos[1]-c[1])^2+(pos[2]-c[2])^2)
    theta::Float64 = atan(pos[2]-c[2], pos[1]-c[1])
    z::Float64 = pos[3]-c[3]

    return [r, theta, z]
end

"""
    checkIndices(stateIndices, stateSize)

Throw error if indices are invalid

# Arguments
- `stateIndices::Vector{Int64}`: Indices of state vector
- `stateSize::Int64`: Number of state get_elements_by_tagname
"""
function checkIndices(stateIndices::Vector{Int64}, stateSize::Int64)
    state::Vector{Float64} = zeros(Float64, stateSize)
    for i::Int64 in 1:length(stateIndices)
        ((stateIndices[i] < 1) || (stateIndices[i] > stateSize)) && throw(BoundsError(state, stateIndices[i]))
        for j in 1:length(stateIndices)
            ((i != j) && (stateIndices[i] == stateIndices[j])) && throw(ArgumentError("Index $(stateIndices[i]) is included more than once"))
        end
    end
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
    n_col::Int64 = length(filter(x -> x == true, mask))
    (size(data, 2) == length(mask)) || throw(ArgumentError("Mask length, $(length(mask)), must match data row length, $(size(data, 2))"))
    if n_col == length(mask)
        return data
    elseif n_col == 0
        return zeros(Float64, (1, 0))
    else
        count::Int64 = 1
        maskedData::Matrix{Float64} = zeros(Float64, (size(data, 1), n_col))
        for c::Int64 in 1:length(mask)
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
