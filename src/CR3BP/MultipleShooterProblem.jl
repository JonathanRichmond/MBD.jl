"""
CR3BP multiple shooter problem wrapper

Author: Jonathan Richmond
C: 9/7/22
U: 1/15/25
"""

import StaticArrays
import MBD: CR3BPMultipleShooterProblem, UNINITIALIZED_INDEX

export addConstraint!, addSegment!, addVariable!, buildAdjacencyMatrix, buildProblem!
export checkJacobian, checkValidGraph, getConstraints, getConstraintVector!
export getFreeVariableIndexMap!, getFreeVariableVector!, getJacobian!, getNumConstraints
export getNumFreeVariables!, importFreeVariables!, removeConstraint!, resetPropagatedArcs!
export setFreeVariableVector!, updateConstraintIndexMap!, updateFreeVariableIndexMap!

"""
    addConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `constraint::AbstractConstraint`: Constraint
"""
function addConstraint!(multipleShooterProblem::CR3BPMultipleShooterProblem, constraint::MBD.AbstractConstraint)
    multipleShooterProblem.constraintIndexMap[constraint] = UNINITIALIZED_INDEX
    updateConstraintIndexMap!(multipleShooterProblem)
end

"""
    addSegment!(multipleShooterProblem, segment)

Return multiple shooter problem object with segment

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `segment::Segment`: Segment
"""
function addSegment!(multipleShooterProblem::CR3BPMultipleShooterProblem, segment::MBD.CR3BPSegment)
    push!(multipleShooterProblem.segments, segment)
    multipleShooterProblem.hasBeenBuilt = false
end

"""
    addVariable!(multipleShooterProblem, variable)

Return multiple shooter problem object with variable

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `variable::Variable`: Variable
"""
function addVariable!(multipleShooterProblem::CR3BPMultipleShooterProblem, variable::MBD.Variable)
    multipleShooterProblem.freeVariableIndexMap[variable] = UNINITIALIZED_INDEX
    updateFreeVariableIndexMap!(multipleShooterProblem)
    resetPropagatedArcs!(multipleShooterProblem)
end

"""
    buildAdjacencyMatrix(multipleShooterProblem)

Return node and segment adjacency matrix

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function buildAdjacencyMatrix(multipleShooterProblem::CR3BPMultipleShooterProblem)
    numNodes::Int16 = length(multipleShooterProblem.nodes)
    adjacencyMatrix::Matrix{Int16} = zeros(Int16, (numNodes,numNodes))
    for segmentIndex::Int16 in Int16(1):Int16(length(multipleShooterProblem.segments))
        segment::MBD.CR3BPSegment = multipleShooterProblem.segments[segmentIndex]
        node0::MBD.CR3BPNode = segment.originNode
        nodef::MBD.CR3BPNode = segment.terminalNode
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function buildProblem!(multipleShooterProblem::CR3BPMultipleShooterProblem)
    empty!(multipleShooterProblem.freeVariableIndexMap)
    empty!(multipleShooterProblem.nodes)
    for segment::MBD.CR3BPSegment in multipleShooterProblem.segments
        node0::MBD.CR3BPNode = segment.originNode
        nodef::MBD.CR3BPNode = segment.terminalNode
        node0Exists::Bool = false
        nodefExists::Bool = false
        for node::MBD.CR3BPNode in multipleShooterProblem.nodes
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
    checkJacobian(multipleShooterProblem)

Return true if Jacobian is accurate

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function checkJacobian(multipleShooterProblem::CR3BPMultipleShooterProblem)
    stepSize::Float64 = sqrt(eps(Float64))
    relTol::Float64 = 2E-3
    problem::CR3BPMultipleShooterProblem = shallowClone(multipleShooterProblem)
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
    constraintIndexMap::Dict{MBD.AbstractConstraint, Int16} = problem.constraintIndexMap
    freeVariableIndexMap::Dict{MBD.Variable, Int16} = getFreeVariableIndexMap!(problem)
    reverseConstraintIndexMap::Dict{Int16, MBD.AbstractConstraint} = Dict{Int16, MBD.AbstractConstraint}()
    reverseFreeVariableIndexMap::Dict{Int16, MBD.Variable} = Dict{Int16, MBD.Variable}()
    for (index::MBD.Variable, value::Int16) in freeVariableIndexMap
        numEntries::Int16 = Int16(getNumFreeVariables(index))
        [reverseFreeVariableIndexMap[value+i] = index for i in 1:numEntries]
    end
    for (index::MBD.AbstractConstraint, value::Int16) in constraintIndexMap
        numEntries::Int16 = Int16(getNumConstraintRows(index))
        [reverseConstraintIndexMap[value+i] = index for i in 1:numEntries]
    end
    println(collect(keys(reverseConstraintIndexMap)))
    absDiff::StaticArrays.SMatrix{numConstraints, numFreeVariables, Float64} = StaticArrays.SMatrix{numConstraints, numFreeVariables, Float64}(jacobianNumerical-jacobianAnalytical)
    relDiff::StaticArrays.MMatrix{numConstraints, numFreeVariables, Float64} = StaticArrays.MMatrix{numConstraints, numFreeVariables, Float64}(copy(absDiff))
    for r::Int16 in Int16(1):Int16(numConstraints)
        for c::Int16 in Int16(1):Int16(numFreeVariables)
            (abs(jacobianNumerical[r,c]) > 1E-12) && (relDiff[r,c] = absDiff[r,c]/jacobianNumerical[r,c])
            if abs(relDiff[r,c]) > relTol
                throw(ErrorException("Jacobian error in entry ($r, $c): Expected = $(jacobianNumerical[r,c]); Actual = $(jacobianAnalytical[r,c]) (Relative error = $(relDiff[r,c])); Constraint: $(typeof(reverseConstraintIndexMap[r])), Free Variable: $(typeof(reverseFreeVariableIndexMap[c]))"))
                return false
            end
        end
    end

    return true
end

"""
    checkValidGraph(multipleShooterProblem, adjacencyMatrix)

Return any graph errors

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `adjacencyMatrix::Matrix{Int16}`: Adjacency matrix
"""
function checkValidGraph(multipleShooterProblem::CR3BPMultipleShooterProblem, adjacencyMatrix::Matrix{Int16})
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
    deepClone(mulitpleShooterProblem)

Return deep copy of multiple shooter problem object

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function deepClone(multipleShooterProblem::CR3BPMultipleShooterProblem)
    object = CR3BPMultipleShooterProblem()
    copiedObjectMap::Dict = Dict()
    object.freeVariableIndexMap = Dict{MBD.Variable, Int16}()
    for (index::MBD.Variable, value::Int16) in multipleShooterProblem.freeVariableIndexMap
        variable::MBD.Variable = MBD.deepClone(index)
        copiedObjectMap[hash(index)] = variable
        object.freeVariableIndexMap[variable] = value
    end
    object.nodes = []
    for node::MBD.CR3BPNode in multipleShooterProblem.nodes
        newNode::MBD.CR3BPNode = MBD.shallowClone(node)
        updatePointers!(newNode, copiedObjectMap)
        copiedObjectMap[hash(node)] = newNode
        push!(object.nodes, newNode)
    end
    object.segments = []
    for segment::MBD.CR3BPSegment in multipleShooterProblem.segments
        newSegment::MBD.CR3BPSegment = MBD.shallowClone(segment)
        updatePointers!(newSegment, copiedObjectMap)
        copiedObjectMap[hash(segment)] = newSegment
        push!(object.segments, newSegment)
    end
    object.constraintIndexMap = Dict{MBD.AbstractConstraint, Int16}()
    for (index::MBD.AbstractConstraint, value::Int16) in multipleShooterProblem.constraintIndexMap
        constraint::MBD.AbstractConstraint = MBD.shallowClone(index)
        updatePointers!(constraint, copiedObjectMap)
        object.constraintIndexMap[constraint] = value
    end
    object.freeVariableVector = copy(multipleShooterProblem.freeVariableVector)
    object.constraintVector = copy(multipleShooterProblem.constraintVector)

    return object
end

"""
    getConstraints(mulitpleShooterProblem)

Return constraints

# Arguments
- `multipleShooterProblem::CR3BPMultipelShooterProblem`: CR3BP multiple shooter problem object
"""
function getConstraints(multipleShooterProblem::CR3BPMultipleShooterProblem)
    return keys(multipleShooterProblem.constraintIndexMap)
end

"""
    getConstraintVector!(multipleShooterProblem)

Return constraint vector

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getConstraintVector!(multipleShooterProblem::CR3BPMultipleShooterProblem)
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getFreeVariableIndexMap!(multipleShooterProblem::CR3BPMultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    
    return multipleShooterProblem.freeVariableIndexMap
end

"""
    getFreeVariableVector!(multipleShooterProblem)

Return free variable vector

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getFreeVariableVector!(multipleShooterProblem::CR3BPMultipleShooterProblem)
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getJacobian!(multipleShooterProblem::CR3BPMultipleShooterProblem)
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getNumConstraints(multipleShooterProblem::CR3BPMultipleShooterProblem)
    numRows::Int64 = 0
    [numRows += getNumConstraintRows(constraint) for constraint in keys(multipleShooterProblem.constraintIndexMap)]

    return numRows
end

"""
    getNumFreeVariables(multipleShooterProblem)

Return number of free variables

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function getNumFreeVariables!(multipleShooterProblem::CR3BPMultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    numRows::Int64 = 0
    [numRows += getNumFreeVariables(variable) for variable in keys(multipleShooterProblem.freeVariableIndexMap)]

    return numRows
end

"""
    importFreeVariables!(multipleShooterProblem, node)

Return multiple shooter problem object with imported node free variables

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `node::CR3BPNode`: CR3BP node object
"""
function importFreeVariables!(multipleShooterProblem::CR3BPMultipleShooterProblem, node::MBD.CR3BPNode)
    map(var -> addVariable!(multipleShooterProblem, var), getVariables(node))
end

"""
    importFreeVariables!(multipleShooterProblem, segment)

Return multiple shooter problem object with imported segment free variables

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `node::CR3BPSegment`: CR3BP segment object
"""
function importFreeVariables!(multipleShooterProblem::CR3BPMultipleShooterProblem, segment::MBD.CR3BPSegment)
    map(var -> addVariable!(multipleShooterProblem, var), getVariables(segment))
end

"""
    removeConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint removed

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `constraint::AbstractConstraint`: Constraint
"""
function removeConstraint!(multipleShooterProblem::CR3BPMultipleShooterProblem, constraint::MBD.AbstractConstraint)
    delete!(multipleShooterProblem.constraintIndexMap, constraint)
    updateConstraintIndexMap!(multipleShooterProblem)
end

"""
    resetPropagatedArcs!(multipleShooterProblem)

Return multiple shooter problem object with empty arcs

# Arguments
- `multipleSHooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function resetPropagatedArcs!(multipleShooterProblem::CR3BPMultipleShooterProblem)
    map(seg -> resetPropagatedArc!(seg), multipleShooterProblem.segments)
end

"""
    setFreeVariableVector!(multipleShooterProblem, freeVariableVector)

Return multiple shooter problem object with updated free variable vector

# Arguments
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function setFreeVariableVector!(multipleShooterProblem::CR3BPMultipleShooterProblem, freeVariableVector::Vector{Float64})
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function shallowClone(multipleShooterProblem::CR3BPMultipleShooterProblem)
    object = CR3BPMultipleShooterProblem()
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function updateConstraintIndexMap!(multipleShooterProblem::CR3BPMultipleShooterProblem)
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
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function updateFreeVariableIndexMap!(multipleShooterProblem::CR3BPMultipleShooterProblem)
    numFreeVariableRows::Int16 = 1
    for variable::MBD.Variable in keys(multipleShooterProblem.freeVariableIndexMap)
        multipleShooterProblem.freeVariableIndexMap[variable] = numFreeVariableRows
        numFreeVariableRows += getNumFreeVariables(variable)
    end
end
