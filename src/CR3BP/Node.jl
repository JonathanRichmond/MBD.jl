"""
CR3BP node wrapper

Author: Jonathan Richmond
C: 9/5/22
U: 1/15/25
"""

import MBD: CR3BPNode

export getVariables

"""
    getVariables(node)

Return variables

# Arguments
- `node::CR3BPNode`: CR3BP node object
"""
function getVariables(node::CR3BPNode)
    return [node.state, node.epoch]
end

"""
    shallowClone(node)

Return copy of node object

# Arguments
- `node::CR3BPNode`: CR3BP node object
"""
function shallowClone(node::CR3BPNode)
    object = CR3BPNode(node.epoch.data[1], node.state.data, node.dynamicsModel)
    object.epoch = node.epoch
    object.state = node.state
    object.dynamicsModel = node.dynamicsModel

    return object
end

"""
    updatePointers!(node, copiedObjectMap)

Update pointers for node object

# Arguments
- `node::CR3BPNode`: CR3BP node object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(node::CR3BPNode, copiedObjectMap::Dict)
    node.state = updatePointer(node.state, copiedObjectMap, true)
    node.epoch = updatePointer(node.epoch, copiedObjectMap, true)
end
