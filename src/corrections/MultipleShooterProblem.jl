"""
Multiple shooter problem wrapper

Author: Jonathan Richmond
C: 9/7/22
U: 8/8/23
"""

import MBD: MultipleShooterProblem

export addConstraint!, addSegment!, addVariable!, buildAdjacencyMatrix
export buildProblem!, checkJacobian, checkValidGraph, getConstraints
export getConstraintVector!, getFreeVariableIndexMap!, getFreeVariableVector!
export getJacobian, getNumberConstraints, getNumberFreeVariables
export importFreeVariables!, removeConstraint!, resetPropagatedArcs!
export setFreeVariableVector!, updateConstraintIndexMap!
export updateFreeVariableIndexMap!

"""
    addConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint

# Arguments
- `multipleShooterProblem:::MultipleShooterProblem`: Multiple shooter problem object
- `constraint::AbstractConstraint`: Constraint
"""
function addConstraint!(multipleShooterProblem::MultipleShooterProblem, constraint::MBD.AbstractConstraint)
    multipleShooterProblem.constraintIndexMap[constraint] = MBD.UNINITIALIZED_INDEX
    (typeof(constraint) <: MBD.IHasVariables) && importFreeVariables!(multipleShooterProblem, constraint)
    updateConstraintIndexMap!(multipleShooterProblem)
end

"""
    addSegment!(multipleShooterProblem, segment)

Return multiple shooter problem object with segment

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `segment::Segment`: Segment
"""
function addSegment!(multipleShooterProblem::MultipleShooterProblem, segment::MBD.Segment)
    push!(multipleShooterProblem.segments, segment)
    multipleShooterProblem.hasBeenBuilt = false
end

"""
    addVariable!(multipleShooterProblem, variable)

Return multiple shooter problem object with variable

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `variable::Variable`: Variable
"""
function addVariable!(multipleShooterProblem::MultipleShooterProblem, variable::MBD.Variable)
    multipleShooterProblem.freeVariableIndexMap[variable] = MBD.UNINITIALIZED_INDEX
    updateFreeVariableIndexMap!(multipleShooterProblem)
    resetPropagatedArcs!(multipleShooterProblem)
end

"""
    buildAdjacencyMatrix(multipleShooterProblem)

Return node and segment adjacency matrix

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function buildAdjacencyMatrix(multipleShooterProblem::MultipleShooterProblem)
    n_nodes::Int64 = length(multipleShooterProblem.nodes)
    adjacencyMatrix::Matrix{Int64} = zeros(Int64, (n_nodes, n_nodes))
    for segmentIndex::Int64 in 1:length(multipleShooterProblem.segments)
        segment::MBD.Segment = multipleShooterProblem.segments[segmentIndex]
        node0::MBD.Node = segment.originNode
        nodef::MBD.Node = segment.terminalNode
        index0::Int64 = 0
        indexf::Int64 = 0
        for nodeIndex::Int64 in 1:length(multipleShooterProblem.nodes)
            (multipleShooterProblem.nodes[nodeIndex] == node0) && (index0 = nodeIndex)
            (multipleShooterProblem.nodes[nodeIndex] == nodef) && (indexf = nodeIndex)
        end
        (index0 == 0) && throw(ErrorException("Could not find node0 in nodes vector"))
        (indexf == 0) && throw(ErrorException("Could not find nodef in nodes vector"))
        adjacencyMatrix[index0, indexf] = segmentIndex
    end

    return adjacencyMatrix
end

"""
    buildProblem!(multipleShooterProblem)

Return built multiple shooter problem object

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function buildProblem!(multipleShooterProblem::MultipleShooterProblem)
    empty!(multipleShooterProblem.freeVariableIndexMap)
    empty!(multipleShooterProblem.nodes)
    for segment::MBD.Segment in multipleShooterProblem.segments
        node0::MBD.Node = segment.originNode
        nodef::MBD.Node = segment.terminalNode
        node0Exists::Bool = false
        nodefExists::Bool = false
        for node::MBD.Node in multipleShooterProblem.nodes
            (hash(node) == hash(node0)) && (node0Exists = true)
            (hash(node) == hash(nodef)) && (nodefExists = true)
        end
        node0Exists || push!(multipleShooterProblem.nodes, node0)
        nodefExists || push!(multipleShooterProblem.nodes, nodef)
        importFreeVariables!(multipleShooterProblem, node0)
        importFreeVariables!(multipleShooterProblem, nodef)
        importFreeVariables!(multipleShooterProblem, segment)
    end
    for constraint::MBD.AbstractConstraint in keys(multipleShooterProblem.constraintIndexMap)
        (typeof(constraint) <: MBD.IHasVariables) && importFreeVar!(multipleShooterProblem, constraint)
    end
    adjacencyMatrix::Matrix{Int64} = buildAdjacencyMatrix(multipleShooterProblem)
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
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function checkJacobian(multipleShooterProblem::MultipleShooterProblem)
    stepSize::Float64 = 1E-8
    relTol::Float64 = 2E-3
    problem::MultipleShooterProblem = shallowClone!(multipleShooterProblem)
    jacobianNumerical::Matrix{Float64} = zeros(Float64, (getNumberConstraints(problem), getNumberFreeVariables(problem)))
    jacobianAnalytical::Matrix{Float64} = copy(getJacobian(problem))
    for varIndex::Int64 in 1:getNumberFreeVariables(problem)
        perturbedFreeVariables::Vector{Float64} = copy(getFreeVariableVector!(problem))
        perturbedFreeVariables[varIndex] -= stepSize
        setFreeVariableVector!(problem, perturbedFreeVariables)
        constraintVectorMinus::Vector{Float64} = copy(getConstraintVector!(problem))
        perturbedFreeVariables[varIndex] += 2*stepSize
        setFreeVariableVector!(problem, perturbedFreeVariables)
        constraintVectorPlus::Vector{Float64} = copy(getConstraintVector!(problem))
        jacobianNumerical[:,varIndex] = (constraintVectorPlus-constraintVectorMinus)./(2*stepSize)
    end
    constraintIndexMap::Dict{MBD.AbstractConstraint, Int64} = problem.constraintIndexMap
    freeVariableIndexMap::Dict{MBD.Variable, Int64} = getFreeVariableIndexMap!(problem)
    reverseConstraintIndexMap::Dict{Int64, MBD.AbstractConstraint} = Dict{Int64, MBD.AbstractConstraint}()
    reverseFreeVariableIndexMap::Dict{Int64, MBD.Variable} = Dict{Int64, MBD.Variable}()
    for (index::MBD.Variable, value::Int64) in freeVariableIndexMap
        n_entries::Int64 = getNumberFreeVariables(index)
        [reverseFreeVariableIndexMap[value+i] = index for i in 1:n_entries]
    end
    for (index::MBD.AbstractConstraint, value::Int64) in constraintIndexMap
        n_entries::Int64 = getNumberConstraintRows(index)
        [reverseConstraintIndexMap[value+i] = index for i in 1:n_entries]
    end
    absDiff::Matrix{Float64} = jacobianNumerical-jacobianAnalytical
    relDiff::Matrix{Float64} = copy(absDiff)
    isEqual::Bool = true
    for r::Int64 in 1:size(relDiff, 1)
        for c::Int64 in 1:size(relDiff, 2)
            (abs(jacobianNumerical[r,c]) > 1E-12) && (relDiff[r,c] = absDiff[r,c]/jacobianNumerical[r,c])
            if abs(relDiff[r,c]) > relTol
                throw(ErrorException("Jacobian error in entry ($r, $c): Expected = $(jacobianNumerical[r,c]); Actual = $(jacobianAnalytical[r,c]) (Relative error = $(relDiff[r,c]))"))
                isEqual = false
            end
        end
    end

    return isEqual
end

"""
    checkValidGraph(multipleShooterProblem, adjacencyMatrix)

Return any graph errors

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `adjacencyMatrix::Matrix{Int64}`: Adjacency matrix
"""
function checkValidGraph(multipleShooterProblem::MultipleShooterProblem, adjacencyMatrix::Matrix{Int64})
    errors::Vector{String} = []
    nodeIsLinked::Vector{Bool} = Vector{Bool}(undef, size(adjacencyMatrix, 1))
    for r::Int64 in 1:size(adjacencyMatrix, 1)
        segmentCount::Int64 = 0
        TOFSegmentSum::Int64 = 0
        for c::Int64 in 1:size(adjacencyMatrix, 2)
            if adjacencyMatrix[r,c] > 0
                segmentCount += 1
                TOF::MBD.Variable = multipleShooterProblem.segments[adjacencyMatrix[r,c]].TOF
                TOFSegmentSum += sign(getData(TOF)[1])
                nodeIsLinked[r] = true
                (r == c) && push!(errors, "Segment $(adjacencyMatrix[r,c]) links node $r to itself")
            end
        end
        if segmentCount > 2
            push!(errors, "Origin node $r has $segmentCount linked segments, but can only have two at most")
        elseif segmentCount == 2
            (TOFSegmentSum == 0) || push!(errors, "origin node $r is linked to two segments with segment TOF = $(TOFSegmentSum/2)")
        end
    end
    for c::Int64 in 1:size(adjacencyMatrix, 2)
        segmentCount::Int64 = 0
        for r::Int64 in 1:size(adjacencyMatrix, 1)
            if adjacencyMatrix[r,c] > 0
                segmentCount += 1
                nodeIsLinked[c] = true
            end
        end
        (segmentCount > 1) && push!(errors, "Terminal node $c has $segmentCount linked segments, but can only have one at most")
    end
    [push!(errors, "Node $r is not linked to any segments") for r in 1:length(nodeIsLinked) if !nodeIsLinked[r]]

    return errors
end

"""
    deepClone(mulitpleShooterProblem)

Return deep copy of multiple shooter problem object

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function deepClone(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    object = MultipleShooterProblem()
    copiedObjectMap::Dict = Dict()
    object.freeVariableIndexMap = Dict{MBD.Variable, Int64}()
    for (index::MBD.Variable, value::Int64) in multipleShooterProblem.freeVariableIndexMap
        variable::MBD.Variable = MBD.deepClone(index)
        copiedObjectMap[hash(index)] = variable
        object.freeVariableIndexMap[variable] = value
    end
    object.nodes = []
    for node::MBD.Node in multipleShooterProblem.nodes
        newNode::MBD.Node = MBD.shallowClone(node)
        updatePointers!(newNode, copiedObjectMap)
        copiedObjectMap[hash(node)] = newNode
        push!(object.nodes, newNode)
    end
    object.segments = []
    for segment::MBD.Segment in multipleShooterProblem.segments
        newSegment::MBD.Segment = MBD.shallowClone(segment)
        updatePointers!(newSegment, copiedObjectMap)
        copiedObjectMap[hash(segment)] = newSegment
        push!(object.segments, newSegment)
    end
    object.constraintIndexMap = Dict{MBD.AbstractConstraint, Int64}()
    for (index::MBD.AbstractConstraint, value::Int64) in multipleShooterProblem.constraintIndexMap
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
- `multipleShooterProblem::MultipelShooterProblem`: Multiple shooter problem object
"""
function getConstraints(multipleShooterProblem::MultipleShooterProblem)
    return keys(multipleShooterProblem.constraintIndexMap)
end

"""
    getConstraintVector!(multipleShooterProblem)

Return constraint vector

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getConstraintVector!(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.constraintVector = Vector{Float64}(undef, getNumberConstraints(multipleShooterProblem))
    for (index::MBD.AbstractConstraint, value::Int64) in multipleShooterProblem.constraintIndexMap
        data::Vector{Float64} = evaluateConstraint(index, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector)
        multipleShooterProblem.constraintVector[value:value+length(data)-1] = data
    end

    return multipleShooterProblem.constraintVector
end

"""
    getFreeVariableIndexMap!(multipleShooterProblem)

Return free variable index map

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getFreeVariableIndexMap!(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    
    return multipleShooterProblem.freeVariableIndexMap
end

"""
    getFreeVariableVector!(multipleShooterProblem)

Return free variable vector

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getFreeVariableVector!(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    multipleShooterProblem.freeVariableVector = Vector{Float64}(undef, getNumberFreeVariables(multipleShooterProblem))
    for (index::MBD.Variable, value::Int64) in multipleShooterProblem.freeVariableIndexMap
        data::Vector{Float64} = getFreeVariableData(index)
        multipleShooterProblem.freeVariableVector[value:value+length(data)-1] = data
    end

    return multipleShooterProblem.freeVariableVector
end

"""
    getJacobian(multipleShooterProblem)

Return Jacobian matrix

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getJacobian(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    jacobian::Matrix{Float64} = zeros(Float64, (getNumberConstraints(multipleShooterProblem), getNumberFreeVariables(multipleShooterProblem)))
    for (index::MBD.AbstractConstraint, value::Int64) in multipleShooterProblem.constraintIndexMap
        partials::Dict{MBD.Variable, Matrix{Float64}} = getPartials_ConstraintWRTVariables(index, multipleShooterProblem.freeVariableIndexMap, multipleShooterProblem.freeVariableVector)
        for (index2::MBD.Variable, value2::Matrix{Float64}) in partials
            maskedData::Matrix{Float64} = maskData(getFreeVariableMask(index2), value2)
            (length(maskedData[1,:]) > 0) && (jacobian[value:value+size(maskedData, 1)-1, multipleShooterProblem.freeVariableIndexMap[index2]:multipleShooterProblem.freeVariableIndexMap[index2]+size(maskedData, 2)-1] = maskedData)
        end
    end

    return jacobian
end

"""
    getNumberConstraints(multipleShooterProblem)

Return number of constraints

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getNumberConstraints(multipleShooterProblem::MultipleShooterProblem)
    n_rows::Int64 = 0
    [n_rows += getNumberConstraintRows(constraint) for constraint in keys(multipleShooterProblem.constraintIndexMap)]

    return n_rows
end

"""
    getNumberFreeVariables(multipleShooterProblem)

Return number of free variables

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function getNumberFreeVariables(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    n_rows::Int64 = 0
    [n_rows += getNumberFreeVariables(variable) for variable in keys(multipleShooterProblem.freeVariableIndexMap)]

    return n_rows
end

"""
    importFreeVariables!(multipleShooterProblem, object)

Return multiple shooter problem object with free variables

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `object::IHasVariables`: Object
"""
function importFreeVariables!(multipleShooterProblem::MultipleShooterProblem, object::MBD.IHasVariables)
    map(var -> addVariable!(multipleShooterProblem, var), getVariables(object))
end

"""
    removeConstraint!(multipleShooterProblem, constraint)

Return multiple shooter problem object with constraint removed

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: MultipleShooterProblem
- `constraint::AbstractConstraint`: Constraint
"""
function removeConstraint!(multipleShooterProblem::MultipleShooterProblem, constraint::MBD.AbstractConstraint)
    delete!(multipleShooterProblem.constraintIndexMap, constraint)
    updateConstraintIndexMap!(multipleShooterProblem.constraintIndexMap)
end

"""
    resetPropagatedArcs!(multipleShooterProblem)

Return multiple shooter problem object with empty arcs

# Arguments
- `multipleSHooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function resetPropagatedArcs!(multipleShooterProblem::MultipleShooterProblem)
    map(seg -> resetPropagatedArc!(seg), multipleShooterProblem.segments)
end

"""
    setFreeVariableVector!(multipleShooterProblem, freeVariableVector)

Return multiple shooter problem object with updated free variable vector

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `freeVariableVector::Vector{Float64}`: Free variable vector
"""
function setFreeVariableVector!(multipleShooterProblem::MultipleShooterProblem, freeVariableVector::Vector{Float64})
    multipleShooterProblem.freeVariableVector = freeVariableVector
    for (index::MBD.Variable, value::Int64) in multipleShooterProblem.freeVariableIndexMap
        n_rows::Int64 = getNumberFreeVariables(index)
        if n_rows > 0
            freeVariables::Vector{Float64} = freeVariableVector[value:value+n_rows-1]
            setFreeVariableData!(index, freeVariables)
        end
    end
    resetPropagatedArcs!(multipleShooterProblem)
end

"""
    shallowClone!(multipleShooterProblem)

Return copy of multiple shooter problem object

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function shallowClone!(multipleShooterProblem::MultipleShooterProblem)
    multipleShooterProblem.hasBeenBuilt || buildProblem!(multipleShooterProblem)
    object = MultipleShooterProblem()
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
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function updateConstraintIndexMap!(multipleShooterProblem::MultipleShooterProblem)
    numConstraintRows::Int64 = 1
    for constraint::MBD.AbstractConstraint in keys(multipleShooterProblem.constraintIndexMap)
        multipleShooterProblem.constraintIndexMap[constraint] = numConstraintRows
        numConstraintRows += getNumberConstraintRows(constraint)
    end
end

"""
    updateFreeVariableIndexMap!(multipleShooterProblem)

Return multiple shooter problem object with updated free variable indices

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function updateFreeVariableIndexMap!(multipleShooterProblem::MultipleShooterProblem)
    numFreeVariableRows::Int64 = 1
    for variable::MBD.Variable in keys(multipleShooterProblem.freeVariableIndexMap)
        multipleShooterProblem.freeVariableIndexMap[variable] = numFreeVariableRows
        numFreeVariableRows += getNumberFreeVariables(variable)
    end
end
