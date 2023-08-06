"""
Utility functions

Author: Jonathan Richmond
C: 9/7/22
U: 8/5/23
"""

export maskData

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
