"""
Node wrapper

Author: Jonathan Richmond
C: 9/5/22
U: 8/6/23
"""

import MBD: Node

export getVariables

"""
    getVariables(node)

Return variables

# Arguments
- `node::Node`: Node object
"""
function getVariables(node::Node)
    return [node.state, node.epoch]
end

"""
    shallowClone(node)

Return copy of node object

# Arguments
- `node::Node`: Node object
"""
function shallowClone(node::Node)
    object = Node(node.epoch.data[1], node.state.data, node.dynamicsModel)
    object.epoch = node.epoch
    object.state = node.state
    object.dynamicsModel = node.dynamicsModel

    return object
end

"""
    updatePointers!(node, copiedObjectMap)

Return node object with updated pointers

# Arguments
- `node::Node`: Node object
- `copiedObjectMap::Dict`: Map between old and new objects
"""
function updatePointers!(node::Node, copiedObjectMap::Dict)
    node.state = updatePointer(node.state, copiedObjectMap, true)
    node.epoch = updatePointer(node.epoch, copiedObjectMap, true)
end
