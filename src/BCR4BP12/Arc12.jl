"""
BCR4BP P1-P2 arc wrapper

Author: Jonathan Richmond
C: 2/19/25
"""

import MBD: BCR4BP12Arc

export deleteStateAndTime!, getStateByIndex, getStateCount, getTimeByIndex

"""
    deleteStateAndTime!(arc, index)

Return arc object with data corresponding to specified index removed

# Arguments
- `arc::BCR4BP12Arc`: BCR4BP P1-P2 arc object
- `index::Int64`: Element index
"""
function deleteStateAndTime!(arc::BCR4BP12Arc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))
    newIndex::Int64 = (index < 0) ? getStateCount(arc) : index
    deleteat!(arc.states, newIndex)
    deleteat!(arc.times, newIndex)
end

"""
    getStateByIndex(arc, index)

Return state at specified index

# Arguments
- `arc::BCR4BP12Arc`: BCR4BP P1-P2 arc object
- `index::Int64`: Element index
"""
function getStateByIndex(arc::BCR4BP12Arc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))

    (index < 0) ? (return copy(arc.states[end+1+index])) : (return copy(arc.states[index]))
end

"""
    getStateCount(arc)

Return number of elements in arc object

# Arguments
- `arc::BCR4BP12Arc`: BCR4BP P1-P2 arc object
"""
function getStateCount(arc::BCR4BP12Arc)
    return length(arc.states)
end

"""
    getTimeByIndex(arc, index)

Return time at specified index

# Arguments
- `arc::BCR4BP12Arc`: BCR4BP P1-P2 arc object
- `index::Int64`: Element index
"""
function getTimeByIndex(arc::BCR4BP12Arc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.times, index))

    (index < 0) ? (return arc.times[end+1+index]) : (return arc.times[index])
end
