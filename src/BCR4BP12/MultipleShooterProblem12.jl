"""
BCR4BP P1-P2 multiple shooter problem wrapper

Author: Jonathan Richmond
C: 4/9/25
U: 4/22/25
"""

import StaticArrays
import MBD: BCR4BP12MultipleShooterProblem, UNINITIALIZED_INDEX

export addConstraint!, addSegment!, addVariable!, buildAdjacencyMatrix, buildProblem!
export checkJacobian, checkValidGraph, getConstraints, getConstraintVector!
export getFreeVariableIndexMap!, getFreeVariableVector!, getJacobian!, getNumConstraints
export getNumFreeVariables!, importFreeVariables!, removeConstraint!, resetPropagatedArcs!
export setFreeVariableVector!, updateConstraintIndexMap!, updateFreeVariableIndexMap!

"""
    addConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `constraint::AbstractConstraint`: Constraint
"""
function addConstraint!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, constraint::MBD.AbstractConstraint)
    multipleShooterProblem.constraintIndexMap[constraint] = UNINITIALIZED_INDEX
    updateConstraintIndexMap!(multipleShooterProblem)
end

"""
    addSegment!(multipleShooterProblem, segment)

Return multiple shooter problem object with segment

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `segment::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function addSegment!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, segment::MBD.BCR4BP12Segment)
    push!(multipleShooterProblem.segments, segment)
    multipleShooterProblem.hasBeenBuilt = false
end

"""
    addVariable!(multipleShooterProblem, variable)

Return multiple shooter problem object with variable

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `variable::Variable`: Variable
"""
function addVariable!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, variable::MBD.Variable)
    multipleShooterProblem.freeVariableIndexMap[variable] = UNINITIALIZED_INDEX
    updateFreeVariableIndexMap!(multipleShooterProblem)
    resetPropagatedArcs!(multipleShooterProblem)
end

"""
    buildAdjacencyMatrix(multipleShooterProblem)

Return node and segment adjacency matrix

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function buildAdjacencyMatrix(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    numNodes::Int16 = length(multipleShooterProblem.nodes)
    adjacencyMatrix::Matrix{Int16} = zeros(Int16, (numNodes,numNodes))
    for segmentIndex::Int16 in Int16(1):Int16(length(multipleShooterProblem.segments))
        segment::MBD.BCR4BP12Segment = multipleShooterProblem.segments[segmentIndex]
        node0::MBD.BCR4BP12Node = segment.originNode
        nodef::MBD.BCR4BP12Node = segment.terminalNode
        index0::Int16 = Int16(0)
        indexf::Int16 = Int16(0)
        for nodeIndex::Int16 in Int16(1):numNodes
            (multipleShooterProblem.nodes[nodeIndex] == node0) && (index0 = nodeIndex)
            (multipleShooterProblem.nodes[nodeIndex] == nodef) && (indexf = nodeIndex)
        end
        (index0 == Int16(0)) && throw(ErrorException("Could not find node0 in nodes vector"))
        (indexf == Int16(0)) && throw(ErrorException("Could not find nodef in nodes vector"))
        adjacencyMatrix[index0, indexf] = segmentIndex
    end

    return adjacencyMatrix
end

"""
    buildProblem!(multipleShooterProblem)

Return built multiple shooter problem object

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function buildProblem!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    empty!(multipleShooterProblem.freeVariableIndexMap)
    empty!(multipleShooterProblem.nodes)
    for segment::MBD.BCR4BP12Segment in multipleShooterProblem.segments
        node0::MBD.BCR4BP12Node = segment.originNode
        nodef::MBD.BCR4BP12Node = segment.terminalNode
        node0Exists::Bool = false
        nodefExists::Bool = false
        for node::MBD.BCR4BP12Node in multipleShooterProblem.nodes
            (hash(node) == hash(node0)) && (node0Exists = true)
            (hash(node) == hash(nodef)) && (nodefExists = true)
        end
        node0Exists || push!(multipleShooterProblem.nodes, node0)
        nodefExists || push!(multipleShooterProblem.nodes, nodef)
        importFreeVariables!(multipleShooterProblem, node0)
        importFreeVariables!(multipleShooterProblem, nodef)
        importFreeVariables!(multipleShooterProblem, segment)
    end
    adjacencyMatrix::Matrix{Int16} = buildAdjacencyMatrix(multipleShooterProblem)
    errors::Vector{String} = checkValidGraph(multipleShooterProblem, adjacencyMatrix)
    if !isempty(errors)
        errorMessage::Vector{String} = ["Directed graph errors: "]
        map(err -> push!(errorMessage, "$err   "), errors)
        push!(errorMessage, "Invalid directed graph for multiple shooting problem")
        throw(ErrorException(errorMessage))
    end
    updateFreeVariableIndexMap!(multipleShooterProblem)
    resetPropagatedArcs!(multipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt = true
end

"""
    checkJacobian(multipleShooterProblem; relTol)

Return true if Jacobian is accurate

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `relTol::Float64`: Relative tolerance (default = 2E-3)
"""
function checkJacobian(multipleShooterProblem::BCR4BP12MultipleShooterProblem, relTol::Float64 = 1E-2)
    stepSize::Float64 = sqrt(eps(Float64))
    problem::BCR4BP12MultipleShooterProblem = shallowClone(multipleShooterProblem)
    numConstraints::Int64 = getNumConstraints(problem)
    numFreeVariables::Int64 = getNumFreeVariables!(problem)
    jacobianNumerical::StaticArrays.MMatrix{numConstraints, numFreeVariables, Float64} = StaticArrays.MMatrix{numConstraints, numFreeVariables, Float64}(zeros(Float64, (numConstraints, numFreeVariables)))
    jacobianAnalytical::StaticArrays.SMatrix{numConstraints, numFreeVariables, Float64} = StaticArrays.SMatrix{numConstraints, numFreeVariables, Float64}(getJacobian!(problem))
    for varIndex::Int16 in Int16(1):numFreeVariables
        perturbedFreeVariables::Vector{Float64} = copy(getFreeVariableVector!(problem))
        perturbedFreeVariables[varIndex] -= stepSize
        setFreeVariableVector!(problem, perturbedFreeVariables)
        constraintVectorMinus::StaticArrays.SVector{numConstraints, Float64} = StaticArrays.SVector{numConstraints, Float64}(copy(getConstraintVector!(problem)))
        perturbedFreeVariables[varIndex] += 2*stepSize
        setFreeVariableVector!(problem, perturbedFreeVariables)
        constraintVectorPlus::StaticArrays.SVector{numConstraints, Float64} = StaticArrays.SVector{numConstraints, Float64}(copy(getConstraintVector!(problem)))
        jacobianNumerical[:,varIndex] = (constraintVectorPlus-constraintVectorMinus)./(2*stepSize)
    end
    freeVariableIndexMap::Dict{MBD.Variable, Int16} = getFreeVariableIndexMap!(problem)
    reverseConstraintIndexMap::Dict{Int16, MBD.AbstractConstraint} = Dict{Int16, MBD.AbstractConstraint}()
    reverseFreeVariableIndexMap::Dict{Int16, MBD.Variable} = Dict{Int16, MBD.Variable}()
    for (index::MBD.Variable, value::Int16) in freeVariableIndexMap
        numEntries::Int16 = Int16(getNumFreeVariables(index))
        [reverseFreeVariableIndexMap[value+i-1] = index for i in 1:numEntries]
    end
    for (index::MBD.AbstractConstraint, value::Int16) in problem.constraintIndexMap
        numEntries::Int16 = Int16(getNumConstraintRows(index))
        [reverseConstraintIndexMap[value+i-1] = index for i in 1:numEntries]
    end
    absDiff::StaticArrays.SMatrix{Int64(numStates), Int64(numStates), Float64} = STMNumerical.-STMAnalytical
    for r::Int16 in Int16(1):numStates, c::Int16 in Int16(1):numStates
        analytical::Float64 = STMAnalytical[r,c]
        numerical::Float64 = STMNumerical[r,c]
        diff::Float64 = absDiff[r,c]
        useAbs::Bool = ((abs(analytical) < stepSize*1E3) || (abs(numerical) < 1E-12))
        relDiff::Float64 = useAbs ? diff : (diff/abs(numerical))
        errorType::String = useAbs ? "Absolute" : "Relative"
        if relDiff > relTol
            throw(ErrorException("Jacobian error in entry ($r, $c): Expected = $numerical; Actual = $analytical; Difference = $diff; Error = $relDiff ($errorType); Constraint (sub-index: $(r-problem.constraintIndexMap[reverseConstraintIndexMap[r]]+1)) = $(typeof(reverseConstraintIndexMap[r])); Free Variable (sub-index: $(c-freeVariableIndexMap[reverseFreeVariableIndexMap[c]]+1)) = $(reverseFreeVariableIndexMap[c].name)"))
        end
    end

    return true
end

"""
    checkValidGraph(multipleShooterProblem, adjacencyMatrix)

Return any graph errors

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `adjacencyMatrix::Matrix{Int16}`: Adjacency matrix
"""
function checkValidGraph(multipleShooterProblem::BCR4BP12MultipleShooterProblem, adjacencyMatrix::Matrix{Int16})
    errors::Vector{String} = []
    s_adjacency::StaticArrays.SVector{2, Int16} = StaticArrays.SVector(Int16(size(adjacencyMatrix, 1)), Int16(size(adjacencyMatrix, 2)))
    nodeIsLinked::Vector{Bool} = Vector{Bool}(undef, s_adjacency[1])
    for r::Int16 in Int16(1):Int16(s_adjacency[1])
        segmentCount::Int16 = 0
        TOFSegmentSum::Int16 = 0
        for c::Int16 in Int16(1):s_adjacency[2]
            if adjacencyMatrix[r,c] > Int16(0)
                segmentCount += 1
                TOF::MBD.Variable = multipleShooterProblem.segments[adjacencyMatrix[r,c]].TOF
                TOFSegmentSum += sign(getData(TOF)[1])
                nodeIsLinked[r] = true
                (r == c) && push!(errors, "Segment $(adjacencyMatrix[r,c]) links node $r to itself")
            end
        end
        if segmentCount > Int16(2)
            push!(errors, "Origin node $r has $segmentCount linked segments, but can only have two at most")
        elseif segmentCount == Int16(2)
            (TOFSegmentSum == Int16(0)) || push!(errors, "origin node $r is linked to two segments with segment TOF = $(TOFSegmentSum/2)")
        end
    end
    for c::Int16 in Int16(1):s_adjacency[2]
        segmentCount::Int16 = 0
        for r::Int16 in Int16(1):s_adjacency[1]
            if adjacencyMatrix[r,c] > Int16(0)
                segmentCount += 1
                nodeIsLinked[c] = true
            end
        end
        (segmentCount > Int16(1)) && push!(errors, "Terminal node $c has $segmentCount linked segments, but can only have one at most")
    end
    [push!(errors, "Node $r is not linked to any segments") for r in 1:length(nodeIsLinked) if !nodeIsLinked[r]]

    return errors
end

"""
    deepClone(multipleShooterProblem)

Return deep copy of multiple shooter problem object

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function deepClone(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    object = BCR4BP12MultipleShooterProblem()
    copiedObjectMap::IdDict{Any, Any} = IdDict{Any, Any}()
    object.freeVariableIndexMap = Dict{MBD.Variable, Int16}()
    for (index::MBD.Variable, value::Int16) in multipleShooterProblem.freeVariableIndexMap
        variable::MBD.Variable = MBD.deepClone(index)
        copiedObjectMap[index] = variable
        object.freeVariableIndexMap[variable] = value
    end
    object.nodes = []
    for node::MBD.BCR4BP12Node in multipleShooterProblem.nodes
        newNode::MBD.BCR4BP12Node = MBD.shallowClone(node)
        updatePointers!(newNode, copiedObjectMap)
        copiedObjectMap[node] = newNode
        push!(object.nodes, newNode)
    end
    object.segments = []
    for segment::MBD.BCR4BP12Segment in multipleShooterProblem.segments
        newSegment::MBD.BCR4BP12Segment = MBD.shallowClone(segment)
        updatePointers!(newSegment, copiedObjectMap)
        copiedObjectMap[segment] = newSegment
        push!(object.segments, newSegment)
    end
    object.constraintIndexMap = Dict{MBD.AbstractConstraint, Int16}()
    for (index::MBD.AbstractConstraint, value::Int16) in multipleShooterProblem.constraintIndexMap
        constraint::MBD.AbstractConstraint = MBD.shallowClone(index, multipleShooterProblem.nodes[1].dynamicsModel)
        updatePointers!(constraint, copiedObjectMap)
        object.constraintIndexMap[constraint] = value
    end
    object.freeVariableVector = copy(multipleShooterProblem.freeVariableVector)
    object.constraintVector = copy(multipleShooterProblem.constraintVector)

    return object
end

"""
    getConstraints(multipleShooterProblem)

Return constraints

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getConstraints(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    return keys(multipleShooterProblem.constraintIndexMap)
end

"""
    getConstraintVector!(multipleShooterProblem)

Return constraint vector

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getConstraintVector!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    multipleShooterProblem.constraintVector = Vector{Float64}(undef, getNumConstraints(multipleShooterProblem))
    for (index::MBD.AbstractConstraint, value::Int16) in multipleShooterProblem.constraintIndexMap
        data::Vector{Float64} = evaluateConstraint(index, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector)
        multipleShooterProblem.constraintVector[value:value+length(data)-1] = data
    end

    return multipleShooterProblem.constraintVector
end

"""
    getFreeVariableIndexMap!(multipleShooterProblem)

Return free variable index map

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getFreeVariableIndexMap!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    
    return multipleShooterProblem.freeVariableIndexMap
end

"""
    getFreeVariableVector!(multipleShooterProblem)

Return free variable vector

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getFreeVariableVector!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    multipleShooterProblem.freeVariableVector = Vector{Float64}(undef, getNumFreeVariables!(multipleShooterProblem))
    for (index::MBD.Variable, value::Int16) in multipleShooterProblem.freeVariableIndexMap
        data::Vector{Float64} = getFreeVariableData(index)
        multipleShooterProblem.freeVariableVector[value:value+length(data)-1] = data
    end

    return multipleShooterProblem.freeVariableVector
end

"""
    getJacobian(multipleShooterProblem)

Return Jacobian matrix

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getJacobian!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    jacobian::Matrix{Float64} = zeros(Float64, (getNumConstraints(multipleShooterProblem),getNumFreeVariables!(multipleShooterProblem)))
    for (index::MBD.AbstractConstraint, value::Int16) in multipleShooterProblem.constraintIndexMap
        partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(index, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector)
        for (index2::MBD.Variable, value2::Matrix{Float64}) in partials
            maskedData::Matrix{Float64} = maskData(getFreeVariableMask(index2), value2)
            (length(maskedData[1,:]) > 0) && (jacobian[value:value+size(maskedData, 1)-1, multipleShooterProblem.freeVariableIndexMap[index2]:multipleShooterProblem.freeVariableIndexMap[index2]+size(maskedData, 2)-1] = maskedData)
        end
    end

    return jacobian
end

"""
    getNumConstraints(multipleShooterProblem)

Return number of constraints

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getNumConstraints(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    numRows::Int64 = 0
    [numRows += getNumConstraintRows(constraint) for constraint in keys(multipleShooterProblem.constraintIndexMap)]

    return numRows
end

"""
    getNumFreeVariables(multipleShooterProblem)

Return number of free variables

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function getNumFreeVariables!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    numRows::Int64 = 0
    [numRows += getNumFreeVariables(variable) for variable in keys(multipleShooterProblem.freeVariableIndexMap)]

    return numRows
end

"""
    importFreeVariables!(multipleShooterProblem, node)

Return multiple shooter problem object with imported node free variables

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `node::BCR4BP12Node`: BCR4BP P1-P2 node object
"""
function importFreeVariables!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, node::MBD.BCR4BP12Node)
    map(var -> addVariable!(multipleShooterProblem, var), getVariables(node))
end

"""
    importFreeVariables!(multipleShooterProblem, segment)

Return multiple shooter problem object with imported segment free variables

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `node::BCR4BP12Segment`: BCR4BP P1-P2 segment object
"""
function importFreeVariables!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, segment::MBD.BCR4BP12Segment)
    map(var -> addVariable!(multipleShooterProblem, var), getVariables(segment))
end

"""
    removeConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint removed

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `constraint::AbstractConstraint`: Constraint
"""
function removeConstraint!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, constraint::MBD.AbstractConstraint)
    delete!(multipleShooterProblem.constraintIndexMap, constraint)
    updateConstraintIndexMap!(multipleShooterProblem)
end

"""
    resetPropagatedArcs!(multipleShooterProblem)

Return multiple shooter problem object with empty arcs

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function resetPropagatedArcs!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    map(seg -> resetPropagatedArc!(seg), multipleShooterProblem.segments)
end

"""
    setFreeVariableVector!(multipleShooterProblem, freeVariableVector)

Return multiple shooter problem object with updated free variable vector

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function setFreeVariableVector!(multipleShooterProblem::BCR4BP12MultipleShooterProblem, freeVariableVector::Vector{Float64})
    multipleShooterProblem.freeVariableVector = freeVariableVector
    for (index::MBD.Variable, value::Int16) in multipleShooterProblem.freeVariableIndexMap
        numRows::Int16 = Int16(getNumFreeVariables(index))
        if numRows > Int16(0)
            freeVariables::Vector{Float64} = freeVariableVector[value:value+numRows-1]
            setFreeVariableData!(index, freeVariables)
        end
    end
    resetPropagatedArcs!(multipleShooterProblem)
end

"""
    shallowClone(multipleShooterProblem)

Return copy of multiple shooter problem object

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function shallowClone(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    object = BCR4BP12MultipleShooterProblem()
    object.constraintIndexMap = copy(multipleShooterProblem.constraintIndexMap)
    object.freeVariableIndexMap = copy(multipleShooterProblem.freeVariableIndexMap)
    object.freeVariableVector = copy(multipleShooterProblem.freeVariableVector)
    object.constraintVector = copy(multipleShooterProblem.constraintVector)
    object.nodes = copy(multipleShooterProblem.nodes)
    object.segments = copy(multipleShooterProblem.segments)

    return object
end

"""
    updateConstraintIndexMap!(multipleShooterProblem)

Return multiple shooter problem object with updated constraint indices

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function updateConstraintIndexMap!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    numConstraintRows::Int16 = 1
    for constraint::MBD.AbstractConstraint in keys(multipleShooterProblem.constraintIndexMap)
        multipleShooterProblem.constraintIndexMap[constraint] = numConstraintRows
        numConstraintRows += getNumConstraintRows(constraint)
    end
end

"""
    updateFreeVariableIndexMap!(multipleShooterProblem)

Return multiple shooter problem object with updated free variable indices

# Arguments
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function updateFreeVariableIndexMap!(multipleShooterProblem::BCR4BP12MultipleShooterProblem)
    numFreeVariableRows::Int16 = 1
    for variable::MBD.Variable in keys(multipleShooterProblem.freeVariableIndexMap)
        multipleShooterProblem.freeVariableIndexMap[variable] = numFreeVariableRows
        numFreeVariableRows += getNumFreeVariables(variable)
    end
end
