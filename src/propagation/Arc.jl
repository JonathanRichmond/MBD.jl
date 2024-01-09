"""
Arc wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 8/5/23
"""

import MBD: Arc

export deleteStateAndTime!, getStateByIndex, getStateCount, getTimeByIndex
export setParameters!

"""
    deleteStateAndTime!(arc, index)

Return arc object with data corresponding to specified index removed

# Arguments
- `arc::Arc`: Arc object
- `index::Int64`: Element index
"""
function deleteStateAndTime!(arc::Arc, index::Int64)
    (index > size(arc.states, 1)) && throw(BoundsError(arc.states, index))
    deleteat!(arc.states, index)
    deleteat!(arc.times, index)
end

"""
    getStateByIndex(arc, index)

Return state at specified index

# Arguments
- `arc::Arc`: Arc object
- `index::Int64`: Element index
"""
function getStateByIndex(arc::Arc, index::Int64)
    (index > size(arc.states, 1)) && throw(BoundsError(arc.states, index))

    (index < 0) ? (return copy(arc.states[end+1+index])) : (return copy(arc.states[index]))
end

"""
    getStateCount(arc)

Return number of elements in arc object

# Arguments
- `arc::Arc`: Arc object
"""
function getStateCount(arc::Arc)
    return size(arc.states, 1)
end

"""
    getTimeByIndex(arc, index)

Return time at specified index

# Arguments
- `arc::Arc`: Arc object
- `index::Int64`: Element index
"""
function getTimeByIndex(arc::Arc, index::Int64)
    (index > size(arc.times, 1)) && throw(BoundsError(arc.times, index))

    (index < 0) ? (return arc.times[end+1+index]) : (return arc.times[index])
end

"""
    setParameters!(arc, params)

Return arc object with propagation parameters

# Arguments
- `arc::Arc`: Arc object
- `params::Vector{Float64}`: System parameters
"""
function setParameters!(arc::Arc, params::Vector{Float64})
    arc.params = copy(params)
end
