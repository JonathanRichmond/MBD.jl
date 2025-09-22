"""
Keplerian arc wrapper

Author: Jonathan Richmond
C: 9/18/25
"""

import MBD: KArc

export deleteStateAndTime!, getMassParameter, getStateByIndex, getStateCount, getTimeByIndex

"""
    deleteStateAndTime!(arc, index)

Return arc object with data corresponding to specified index removed

# Arguments
- `arc::KArc`: Keplerian arc object
- `index::Int64`: Element index
"""
function deleteStateAndTime!(arc::KArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))
    newIndex::Int64 = (index < 0) ? getStateCount(arc) : index
    deleteat!(arc.states, newIndex)
    deleteat!(arc.times, newIndex)
end

"""
    getMassParameter(arc)

Return Keplerian system mass parameter

# Arguments
- `arc::KArc`: Keplerian arc object
"""
function getMassParameter(arc::KArc)
    return getMassRatio(arc.dynamicsModel)
end

"""
    getStateByIndex(arc, index)

Return state at specified index

# Arguments
- `arc::KArc`: Keplerian arc object
- `index::Int64`: Element index
"""
function getStateByIndex(arc::KArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))

    (index < 0) ? (return copy(arc.states[end+1+index])) : (return copy(arc.states[index]))
end

"""
    getStateCount(arc)

Return number of elements in arc object

# Arguments
- `arc::KArc`: Keplerian arc object
"""
function getStateCount(arc::KArc)
    return length(arc.states)
end

"""
    getTimeByIndex(arc, index)

Return time at specified index

# Arguments
- `arc::KArc`: Keplerian arc object
- `index::Int64`: Element index
"""
function getTimeByIndex(arc::KArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.times, index))

    (index < 0) ? (return arc.times[end+1+index]) : (return arc.times[index])
end
