"""
CR3BP arc wrapper

Author: Jonathan Richmond
C: 9/2/22
U: 1/15/25
"""

import MBD: CR3BPArc

export deleteStateAndTime!, getMassRatio, getStateByIndex, getStateCount, getTimeByIndex

"""
    deleteStateAndTime!(arc, index)

Return arc object with data corresponding to specified index removed

# Arguments
- `arc::CR3BPArc`: CR3BP arc object
- `index::Int64`: Element index
"""
function deleteStateAndTime!(arc::CR3BPArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))
    deleteat!(arc.states, index)
    deleteat!(arc.times, index)
end

"""
    getMassRatio(arc)

Return CR3BP system mass ratio

# Arguments
- `arc::CR3BPArc`: CR3BP arc object
"""
function getMassRatio(arc::CR3BPArc)
    return getMassRatio(arc.dynamicsModel)
end

"""
    getStateByIndex(arc, index)

Return state at specified index

# Arguments
- `arc::CR3BPArc`: CR3BP arc object
- `index::Int64`: Element index
"""
function getStateByIndex(arc::CR3BPArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.states, index))

    (index < 0) ? (return copy(arc.states[end+1+index])) : (return copy(arc.states[index]))
end

"""
    getStateCount(arc)

Return number of elements in arc object

# Arguments
- `arc::CR3BPArc`: CR3BP arc object
"""
function getStateCount(arc::CR3BPArc)
    return length(arc.states)
end

"""
    getTimeByIndex(arc, index)

Return time at specified index

# Arguments
- `arc::CR3BPArc`: CR3BP arc object
- `index::Int64`: Element index
"""
function getTimeByIndex(arc::CR3BPArc, index::Int64)
    (index > getStateCount(arc)) && throw(BoundsError(arc.times, index))

    (index < 0) ? (return arc.times[end+1+index]) : (return arc.times[index])
end
