"""
BCR4BP P1-P2 node wrapper

Author: Jonathan Richmond
C: 4/9/25
"""

import MBD: BCR4BP12Node

export getVariables

"""
    getVariables(node)

Return variables

# Arguments
- `node::BCR4BP12Node`: BCR4BP P1-P2 node object
"""
function getVariables(node::BCR4BP12Node)
    return [node.state, node.epoch]
end

"""
    shallowClone(node)

Return copy of node object

# Arguments
- `node::BCR4BP12Node`: BCR4BP P1-P2 node object
"""
function shallowClone(node::BCR4BP12Node)
    object = BCR4BP12Node(node.epoch.data[1], node.state.data, node.dynamicsModel)
    object.epoch = node.epoch
    object.state = node.state
    object.dynamicsModel = node.dynamicsModel

    return object
end

"""
    updatePointers!(node, copiedObjectMap)

Update pointers for node object

# Arguments
- `node::BCR4BP12Node`: BCR4BP P1-P2 node object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(node::BCR4BP12Node, copiedObjectMap::Dict)
    node.state = updatePointer(node.state, copiedObjectMap, true)
    node.epoch = updatePointer(node.epoch, copiedObjectMap, true)
end
