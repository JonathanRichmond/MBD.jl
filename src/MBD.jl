"""
Multi-body dynamics astrodynamics package

Author: Jonathan Richmond
C: 9/1/22
U: 2/19/25
"""
module MBD

import Base: ==
import Combinatorics, DifferentialEquations, LightXML, LinearAlgebra, SPICE, StaticArrays

const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0

"""
Enumerated type for bifurcations
"""
@enum BifurcationType TANGENT DOUBLING TRIPLING QUADRUPLING HOPF

"""
Enumerated type for EOMs
"""
@enum EquationType ARCLENGTH FULL MOMENTUM SIMPLE STM

"""
Enumerated type for integrators
"""
@enum IntegratorType AB5 ABM54 BS5 DP5 DP8

"""
Abstract type for constraints
"""
abstract type AbstractConstraint end

"""
Abstract type for continuation end checks
"""
abstract type AbstractContinuationEndCheck end

"""
Abstract type for continuation jump checks
"""
abstract type AbstractContinuationJumpCheck end

"""
Abstract type for targeters
"""
abstract type AbstractTargeter end

"""
Abstract type for update generators
"""
abstract type AbstractUpdateGenerator end

"""
    BodyName(name)

Body name object

# Arguments
- `name::String`: Name of body
"""
struct BodyName
    name::String                                                        # Name

    function BodyName(name::String)
        this = new(name)

        return this
    end
end

"""
    BodyData(name)

Body data object

# Arguments
- `name::String`: Name of body
"""
mutable struct BodyData
    bodyRadius::Float64                                                 # Body mean radius [km]
    gravParam::Float64                                                  # Gravitational parameter [kg^3/s^2]
    inc::Float64                                                        # Orbit inclination relative to Ecliptic J2000 frame [rad]
    mass::Float64                                                       # Mass [kg]
    name::String                                                        # Name
    orbitRadius::Float64                                                # Orbit mean radius [km]
    parentSpiceID::Int16                                                # Parent body SPICE ID
    RAAN::Float64                                                       # Orbit Right-Ascension of Ascending Node relative to Ecliptic J2000 frame [rad]
    spiceID::Int16                                                      # SPICE ID

    function BodyData(name::String)
        this = new()

        dataFile::String = joinpath(@__DIR__, "body_data.xml")
        bodyName = BodyName(name)
        id::Int16 = getIDCode(bodyName)
        try
            doc::LightXML.XMLDocument = LightXML.parse_file(dataFile)
            root::LightXML.XMLElement = LightXML.root(doc)
            nodeList::Vector{LightXML.XMLElement} = LightXML.get_elements_by_tagname(root, "body")
            for node::LightXML.XMLElement in nodeList
                if LightXML.is_elementnode(node)
                    elID::String = LightXML.content(LightXML.get_elements_by_tagname(node, "id")[1])
                    this.spiceID = parse(Int16, elID)
                    if this.spiceID == id
                        this.name = LightXML.content(LightXML.get_elements_by_tagname(node, "name")[1])
                        parentIDStr::String = LightXML.content(LightXML.get_elements_by_tagname(node, "parentId")[1])
                        (parentIDStr == "NaN") || (this.parentSpiceID = parse(Int16, parentIDStr))
                        this.gravParam = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "gm")[1]))
                        this.mass = this.gravParam/GRAVITY
                        this.bodyRadius = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "radius")[1]))
                        this.orbitRadius = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "circ_r")[1]))
                        this.inc = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "inc")[1]))
                        this.RAAN = parse(Float64, LightXML.content(LightXML.get_elements_by_tagname(node, "raan")[1]))
                        break
                    end
                end
            end
            isdefined(this, :name) || throw(ArgumentError("Body $id not found in 'body_data.xml'"))
        catch
            throw(ErrorException("Could not parse $dataFile"))
        end

        return this
    end
end
Base.:(==)(bodyData1::BodyData, bodyData2::BodyData) = (bodyData1.spiceID == bodyData2.spiceID)

"""
    CR3BPSystemData(p1, p2)

CR3BP system data object

# Arguments
- `p1::String`: Name of first primary
- `p2::String`: Name of second primary
"""
mutable struct CR3BPSystemData
    primaryData::StaticArrays.SVector{2, BodyData}                      # Primary data objects
    primaryNames::StaticArrays.SVector{2, String}                       # Primary names
    primarySpiceIDs::StaticArrays.SVector{2, Int16}                     # Primary SPICE IDs

    function CR3BPSystemData(p1::String, p2::String)
        this = new()
        
        this.primaryData = StaticArrays.SVector(BodyData(p1), BodyData(p2))
        this.primaryNames = StaticArrays.SVector(p1, p2)
        this.primarySpiceIDs = StaticArrays.SVector(this.primaryData[1].spiceID, this.primaryData[2].spiceID)
        (this.primaryData[2].parentSpiceID == this.primarySpiceIDs[1]) || throw(ArgumentError("First primary must be parent of second primary"))

        return this
    end
end
Base.:(==)(systemData1::CR3BPSystemData, systemData2::CR3BPSystemData) = ((systemData1.primaryData == systemData2.primaryData) && (systemData1.primaryNames == systemData2.primaryNames) && (systemData1.primarySpiceIDs == systemData2.primarySpiceIDs))

"""
    CR3BPDynamicsModel(systemData)

CR3BP dynamics model object

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system data object
"""
struct CR3BPDynamicsModel
    systemData::CR3BPSystemData                                         # CR3BP system data object

    function CR3BPDynamicsModel(systemData::CR3BPSystemData)
        this = new(systemData)

        return this
    end
end
Base.:(==)(dynamicsModel1::CR3BPDynamicsModel, dynamicsModel2::CR3BPDynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

"""
    CR3BPEquationsOfMotion(equationType, dynamicsModel)

CR3BP EOM object

# Arguments
- `equationType::EquationType`: EOM type
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct CR3BPEquationsOfMotion
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    equationType::EquationType                                          # EOM type
    
    function CR3BPEquationsOfMotion(equationType::EquationType, dynamicsModel::CR3BPDynamicsModel)
        this = new(dynamicsModel, equationType)

        return this
    end
end
Base.:(==)(EOMs1::CR3BPEquationsOfMotion, EOMs2::CR3BPEquationsOfMotion) = ((EOMs1.dynamicsModel == EOMs2.dynamicsModel) && (EOMs1.equationType == EOMs2.equationType))

"""
    IntegratorFactory(integratorSolver)

Integrator factory object

# Arguments
- `integratorSolver::IntegratorType`: Integrator type
"""
mutable struct IntegratorFactory
    integrator                                                          # Integrator object
    integratorType::IntegratorType                                      # Integrator type
    
    function IntegratorFactory(integratorSolver::IntegratorType)
        this = new()

        this.integratorType = integratorSolver
        if integratorSolver == DP8
            this.integrator = DifferentialEquations.DP8()
        elseif integratorSolver == AB5
            this.integrator = DifferentialEquations.AB5()
        elseif integratorSolver == ABM54
            this.integrator = DifferentialEquations.ABM54()
        elseif integratorSolver == BS5
            this.integrator = DifferentialEquations.BS5()
        elseif integratorSolver == DP5
            this.integrator = DifferentialEquations.DP5()
        else
            throw(ArgumentError("Invalid integrator type"))
        end

        return this
    end
end
Base.:(==)(integratorFactory1::IntegratorFactory, integratorFactory2::IntegratorFactory) = (integratorFactory1.integratorType == integratorFactory2.integratorType)

"""
    Propagator(; integratorSolver, equationType)

Propagator object

# Arguments
- `integratorSolver::IntegratorType`: Integrator type (default = DP8)
- `equationType::EquationType`: Equation type (default = SIMPLE)
"""
mutable struct Propagator
    absTol::Float64                                                     # Absolute tolerance
    equationType::EquationType                                          # EOM type
    integratorFactory::IntegratorFactory                                # Integrator factory object
    maxEvaluationCount::Int64                                           # Maximum number of equation evaluations
    maxStep::Int64                                                      # Maximum step size
    relTol::Float64                                                     # Relative tolerance

    function Propagator(; integratorSolver::IntegratorType = DP8, equationType::EquationType = SIMPLE)
        this = new()

        this.integratorFactory = IntegratorFactory(integratorSolver)
        this.equationType = equationType
        this.absTol = 1E-12
        this.relTol = 1E-12
        this.maxEvaluationCount = typemax(Int64)
        this.maxStep = 100

        return this
    end
end
Base.:(==)(propagator1::Propagator, propagator2::Propagator) = ((propagator1.equationType == propagator2.equationType) && (propagator1.integratorFactory == propagator2.integratorFactory))

"""
    CR3BPArc(dynamicsModel)

CR3BP arc object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
mutable struct CR3BPArc
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    states::Vector{Vector{Float64}}                                     # State vectors along arc [ndim]
    times::Vector{Float64}                                              # Times along arc [ndim]

    function CR3BPArc(dynamicsModel::CR3BPDynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
        this.states = []
        this.times = []

        return this
    end
end
Base.:(==)(arc1::CR3BPArc, arc2::CR3BPArc) = ((arc1.dynamicsModel == arc2.dynamicsModel) && (arc1.states == arc2.states) && (arc1.times == arc2.times))

"""
    Variable(data, freeVarMask)

Variable object

# Arguments
- `data::Vector{Float64}`: Data values
- `freeVarMask::Vector{Bool}`: Free variable mask
"""
mutable struct Variable
    data::Vector{Float64}                                               # Data values
    freeVariableMask::Vector{Bool}                                      # Free variable mask
    name::String                                                        # Name

    function Variable(data::Vector{Float64}, freeVariableMask::Vector{Bool})
        this = new()

        (length(freeVariableMask) == length(data)) || throw(ArgumentError("Free variable mask length, $(length(freeVariableMask)), must match data values length, $(length(data))"))
        this.name = ""
        this.data = copy(data)
        this.freeVariableMask = freeVariableMask

        return this
    end
end
Base.:(==)(variable1::Variable, variable2::Variable) = ((variable1.data == variable2.data) && (variable1.freeVariableMask == variable2.freeVariableMask))

"""
    CR3BPNode(epoch, state, dynamicsModel)

CR3BP node object

# Arguments
- `epoch::Float64`: Epoch [ndim]
- `state::Vector{Float64}`: State vector [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
mutable struct CR3BPNode
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    epoch::Variable                                                     # Epoch variable
    state::Variable                                                     # State variable

    function CR3BPNode(epoch::Float64, state::Vector{Float64}, dynamicsModel::CR3BPDynamicsModel)
        this = new()

        if isEpochIndependent(dynamicsModel)
            this.epoch = Variable([epoch], [false])
        else
            this.epoch = Variable([epoch], [true])
        end
        this.epoch.name = "Node Epoch"
        n_simple::Int64 = getStateSize(dynamicsModel, SIMPLE)
        (length(state) < n_simple) && throw(ArgumentError("State vector length is $(length(state)), but should be $n_simple"))
        this.state = Variable(((length(state) > n_simple) ? copy(state[1:n_simple]) : copy(state)), [true for n in 1:length(state)])
        this.state.name = "Node State"
        this.dynamicsModel = dynamicsModel

        return this
    end
end
Base.:(==)(node1::CR3BPNode, node2::CR3BPNode) = ((node1.dynamicsModel == node2.dynamicsModel) && (node1.epoch == node2.epoch) && (node1.state == node2.state))

"""
    CR3BPSegment(TOF, originNode, terminalNode)

CR3BP segment object

# Arguments
- `TOF::Float64`: Time-of-flight
- `originNode::CR3BPNode`: Origin node
- `terminalNode::CR3BPNode`: Terminal node
"""
mutable struct CR3BPSegment
    originNode::CR3BPNode                                               # Origin CR3BP node object
    propArc::CR3BPArc                                                   # Propagation CR3BP arc object
    propagator::Propagator                                              # Propagator object
    terminalNode::CR3BPNode                                             # Terminal CR3BP node object
    TOF::Variable                                                       # Time-of-flight variable

    function CR3BPSegment(TOF::Float64, originNode::CR3BPNode, terminalNode::CR3BPNode)
        this = new()

        (originNode == terminalNode) && throw(ArgumentError("Origin and terminal nodes cannot be identical"))
        this.TOF = Variable([TOF], [true])
        this.TOF.name = "Segment Time-of-Flight"
        this.originNode = originNode
        this.propArc = CR3BPArc(originNode.dynamicsModel)
        this.terminalNode = terminalNode
        this.propagator = Propagator()

        return this
    end
end
Base.:(==)(segment1::CR3BPSegment, segment2::CR3BPSegment) = ((segment1.originNode == segment2.originNode) && (segment1.terminalNode == segment2.terminalNode) && (segment1.TOF == segment2.TOF))

"""
    CR3BPMultipleShooterProblem()

CR3BP multiple shooter problem object
"""
mutable struct CR3BPMultipleShooterProblem
    constraintIndexMap::Dict{AbstractConstraint, Int16}                 # Map between constraints and first index of contraint equation in constraint vector
    constraintVector::Vector{Float64}                                   # Constraints
    freeVariableIndexMap::Dict{Variable, Int16}                         # Map between free variables and first index of free variables in free variable vector
    freeVariableVector::Vector{Float64}                                 # Free variables
    hasBeenBuilt::Bool                                                  # Has been built?
    jacobian::Vector{Vector{Float64}}                                   # Jacobian matrix
    nodes::Vector{CR3BPNode}                                            # CR3BP node objects
    segments::Vector{CR3BPSegment}                                      # CR3BP segment objects

    function CR3BPMultipleShooterProblem()
        this = new()

        this.nodes = []
        this.segments = []
        this.freeVariableVector = []
        this.freeVariableIndexMap = Dict{Variable, Int16}()
        this.constraintVector = []
        this.constraintIndexMap = Dict{AbstractConstraint, Int16}()
        this.jacobian = []
        this.hasBeenBuilt = false

        return this
    end
end
Base.:(==)(problem1::CR3BPMultipleShooterProblem, problem2::CR3BPMultipleShooterProblem) = ((problem1.constraintVector == problem2.constraintVector) && (problem1.freeVariableVector == problem2.freeVariableVector)  && (problem1.jacobian == problem2.jacobian) && (problem1.nodes == problem2.nodes) && (problem1.segments == problem2.segments))

"""
    CR3BPContinuityConstraint(segment)

CR3BP continuity constraint object

# Arguments
- `segment::CR3BPSegment`: CR3BP segment object
"""
mutable struct CR3BPContinuityConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int16}                                   # Constrained state indices
    segment::CR3BPSegment                                               # CR3BP segment object

    function CR3BPContinuityConstraint(segment::CR3BPSegment)
        this = new()

        this.segment = segment
        this.constrainedIndices = Int16(1):Int16(length(segment.originNode.state.data))

        return this
    end
end
Base.:(==)(continuityConstraint1::CR3BPContinuityConstraint, continuityConstraint2::CR3BPContinuityConstraint) = ((continuityConstraint1.constrainedIndices == continuityConstraint2.constrainedIndices) && (continuityConstraint1.segment == continuityConstraint2.segment))

"""
    CR3BPStateConstraint(node, indices, values)

CR3BP state constraint object

# Arguments
- `node::CR3BPNode`: CR3BP node object
- `indices::Vector{Int64}`: Constrained state indices
- `values::Vector{Float64}`: Constraint values
"""
mutable struct CR3BPStateConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int16}                                   # Constrained state indices
    values::Vector{Float64}                                             # Constraint values
    variable::Variable                                                  # Constrained variable

    function CR3BPStateConstraint(node::CR3BPNode, indices::Vector{Int64}, values::Vector{Float64})
        this = new()

        (length(indices) == length(values)) || throw(ArgumentError("Number of indices, $(length(indices)), must match number of values, $(length(values))"))
        this.variable = node.state
        checkIndices(indices, length(this.variable.data))
        this.constrainedIndices = convert(Vector{Int16}, indices)
        this.values = copy(values)
        freeVariableMask::Vector{Bool} = getFreeVariableMask(this.variable)
        for index::Int64 in this.constrainedIndices
            freeVariableMask[index] || throw(ArgumentError("Variable element $index is constrained but not free to vary"))
        end

        return this
    end
end
Base.:(==)(stateConstraint1::CR3BPStateConstraint, stateConstraint2::CR3BPStateConstraint) = ((stateConstraint1.constrainedIndices == stateConstraint2.constrainedIndices) && (stateConstraint1.values == stateConstraint2.values) && (stateConstraint1.variable == stateConstraint2.variable))

"""
    StateMatchConstraint(state1, state2, indices)

State match constraint object

# Arguments
- `state1::Variable`: First state vector
- `state2::Variable`: Second state vector
- `indices::Vector{Int64}`: Constrained state indices
"""
mutable struct StateMatchConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int16}                                   # Constrained state indices
    variable1::Variable                                                 # First state variable
    variable2::Variable                                                 # Second state variable

    function StateMatchConstraint(state1::Variable, state2::Variable, indices::Vector{Int64})
        this = new()

        (length(state1.data) == length(state2.data)) || throw(ArgumentError("First state length, $(length(state1.data)), must match second state length, $(length(state2.data))"))
        this.variable1 = state1
        this.variable2 = state2
        checkIndices(indices, length(state1.data))
        this.constrainedIndices = convert(Vector{Int16}, indices)
        mask1::Vector{Bool} = getFreeVariableMask(this.variable1)
        mask2::Vector{Bool} = getFreeVariableMask(this.variable2)
        for i::Int64 in eachindex(this.constrainedIndices)
            (mask1[this.constrainedIndices[i]] || mask2[this.constrainedIndices[i]]) || throw(ArgumentError("Cannot constrain index = $(this.constrainedIndices[i]); it is not a free variable"))
        end

        return this
    end
end
Base.:(==)(stateMatchConstraint1::StateMatchConstraint, stateMatchConstraint2::StateMatchConstraint) = ((stateMatchConstraint1.constrainedIndices == stateMatchConstraint2.constrainedIndices) && (stateMatchConstraint1.variable1 == stateMatchConstraint2.variable1) && (stateMatchConstraint1.variable2 == stateMatchConstraint2.variable2))

"""
    JacobiConstraint(node, value)

Jacobi constraint object

# Arguments
- `node::CR3BPNode`: CR3BP node object
- `value::Float64`: Constraint value
"""
mutable struct JacobiConstraint <: AbstractConstraint
    dynamicsModel::CR3BPDynamicsModel                                   # Dynamics model object
    epoch::Variable                                                     # Epoch
    state::Variable                                                     # Constrained state
    value::Float64                                                      # Constraint value

    function JacobiConstraint(node::CR3BPNode, value::Float64)
        this = new()

        (typeof(node.dynamicsModel) == CR3BPDynamicsModel) ? (this.dynamicsModel = node.dynamicsModel) : throw(ArgumentError("Expecting CR3BP dynamics model object"))
        this.epoch = node.epoch
        this.state = node.state
        this.value = value

        return this
    end
end
Base.:(==)(jacobiConstraint1::JacobiConstraint, jacobiConstraint2::JacobiConstraint) = ((jacobiConstraint1.dynamicsModel == jacobiConstraint2.dynamicsModel) && (jacobiConstraint1.epoch == jacobiConstraint2.epoch) && (jacobiConstraint1.state == jacobiConstraint2.state) && (jacobiConstraint1.value == jacobiConstraint2.value))

"""
    ConstraintVectorL2NormConvergenceCheck(tol)

Constraint vector L2 norm convergence check object

# Arguments
- `tol::Float64`: Convergence tolerance (default = 1E-10)
"""
struct ConstraintVectorL2NormConvergenceCheck
    maxVectorNorm::Float64                                              # Maximum allowable vector norm

    function ConstraintVectorL2NormConvergenceCheck(tol::Float64 = 1E-10)
        this = new(tol)

        return this
    end
end
Base.:(==)(constraintVectorL2NormConvergenceCheck1::ConstraintVectorL2NormConvergenceCheck, constraintVectorL2NormConvergenceCheck2::ConstraintVectorL2NormConvergenceCheck) = (constraintVectorL2NormConvergenceCheck1.maxVectorNorm == constraintVectorL2NormConvergenceCheck2.maxVectorNorm)

"""
    MinimumNormUpdateGenerator()

Minimum norm update generator object
"""
struct MinimumNormUpdateGenerator <: AbstractUpdateGenerator
    function MinimumNormUpdateGenerator()
        this = new()

        return this
    end
end

"""
    LeastSquaresUpdateGenerator()

Least squares update generator object
"""
struct LeastSquaresUpdateGenerator <: AbstractUpdateGenerator
    function LeastSquaresUpdateGenerator()
        this = new()

        return this
    end
end

"""
    CR3BPMultipleShooter(tol)

CR3BP multiple shooter object

# Arguments
- `tol::Float64`: Convergence tolerance (default = 1E-10)
"""
mutable struct CR3BPMultipleShooter
    convergenceCheck::ConstraintVectorL2NormConvergenceCheck            # Convergence check object
    maxIterations::Int16                                                # Maximum number of solver iterations
    printProgress::Bool                                                 # Print progress?
    recentIterationCount::Int16                                         # Number of iterations during last solve
    solutionInProgress::CR3BPMultipleShooterProblem                     # CR3BP multiple shooter problem object being solved
    updateGenerators::StaticArrays.SVector{2, AbstractUpdateGenerator}  # Update generator objects

    function CR3BPMultipleShooter(tol::Float64 = 1E-10)
        this = new()

        this.recentIterationCount = Int16(0)
        this.solutionInProgress = CR3BPMultipleShooterProblem()
        this.convergenceCheck = ConstraintVectorL2NormConvergenceCheck(tol)
        this.updateGenerators = StaticArrays.SVector(MinimumNormUpdateGenerator(), LeastSquaresUpdateGenerator())
        this.maxIterations = Int16(25)
        this.printProgress = false


        return this
    end
end
Base.:(==)(multipleShooter1::CR3BPMultipleShooter, multipleShooter2::CR3BPMultipleShooter) = ((multipleShooter1.convergenceCheck == multipleShooter2.convergenceCheck) && (multipleShooter1.maxIterations == multipleShooter2.maxIterations) && (multipleShooter1.recentIterationCount == multipleShooter2.recentIterationCount) && (multipleShooter1.solutionInProgress == multipleShooter2.solutionInProgress))

"""
    CR3BPContinuationFamily(dynamicsModel)

CR3BP continuation family object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
mutable struct CR3BPContinuationFamily
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    nodes::Vector{Vector{CR3BPNode}}                                    # CR3BP node objects
    segments::Vector{Vector{CR3BPSegment}}                              # CR3BP segment objects

    function CR3BPContinuationFamily(dynamicsModel::CR3BPDynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
        this.nodes = []
        this.segments = []

        return this
    end
end
Base.:(==)(continuationFamily1::CR3BPContinuationFamily, continuationFamily2::CR3BPContinuationFamily) = ((continuationFamily1.dynamicsModel == continuationFamily2.dynamicsModel) && (continuationFamily1.nodes == continuationFamily2.nodes) && (continuationFamily1.segments == continuationFamily2.segments))

"""
    CR3BPContinuationData(solution1, solution2)

CR3BP continuation data object

# Arguments
- `solution1::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
- `solution2::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
"""
mutable struct CR3BPContinuationData
    converging::Bool                                                    # Converging?
    currentStepSize::Float64                                            # Current continuation step size
    family::CR3BPContinuationFamily                                     # Computed CR3BP continuation family object
    forceEndContinuation::Bool                                          # Force end of continuation?
    fullStep::Vector{Float64}                                           # Full step along free variable vector
    initialGuess::CR3BPMultipleShooterProblem                           # First member passed to corrections algorithm
    nextGuess::CR3BPMultipleShooterProblem                              # Inital guess for next family member
    numIterations::Int16                                                # Number of iterations required to converge previous solution
    previousSolution::CR3BPMultipleShooterProblem                       # Most recently converged family member
    twoPreviousSolution::CR3BPMultipleShooterProblem                    # Second previously converged family member

    function CR3BPContinuationData(solution1::CR3BPMultipleShooterProblem, solution2::CR3BPMultipleShooterProblem)
        this = new()

        this.previousSolution = solution1
        this.twoPreviousSolution = solution2
        this.numIterations = Int16(0)
        this.initialGuess = CR3BPMultipleShooterProblem()
        this.converging = true
        this.fullStep = []
        this.currentStepSize = 1.0
        this.nextGuess = CR3BPMultipleShooterProblem()
        this.family = CR3BPContinuationFamily(solution1.nodes[1].dynamicsModel)
        this.forceEndContinuation = false

        return this
    end
end
Base.:(==)(continuationData1::CR3BPContinuationData, continuationData2::CR3BPContinuationData) = ((continuationData1.currentStepSize = continuationData2.currentStepSize) && (continuationData1.family = continuationData2.family) && (continuationData1.fullStep = continuationData2.fullStep) && (continuationData1.initialGuess = continuationData2.initialGuess) && (continuationData1.nextGuess = continuationData2.nextGuess) && (continuationData1.numIterations = continuationData2.numIterations) && (continuationData1.previousSolution = continuationData2.previousSolution) && (continuationData1.twoPreviousSolution = continuationData2.twoPreviousSolution))

"""
    AdaptiveStepSizeByElementGenerator(elementName, elementIndex, initialStepSize, maxElementStepSize)

Adaptive step size by element generator object

# Arguments
- `elementName::String`: Free variable name
- `elementIndex::Int64`: Free variable index
- `initialStepSize::Float64`: Initial step size
- `maxElementStepSize::Float64`: Maximum element step size
"""
mutable struct AdaptiveStepSizeByElementGenerator
    elementIndex::Int16                                                 # Free variable index
    elementName::String                                                 # Free variable name
    initialStepSize::Float64                                            # Initial continuation step size
    maxElementStepSize::Float64                                         # Maximum continuation step size for element
    maxIterations::Int16                                                # Maximum iterations to converge solution
    maxStepSize::Float64                                                # Maximum continuation step size
    minIterations::Int16                                                # Minimum iterations to converge solution
    minStepSize::Float64                                                # Minimum continuation step size
    scaleFactor::Float64                                                # Step size scale scaleFactor

    function AdaptiveStepSizeByElementGenerator(elementName::String, elementIndex::Int64, initialStepSize::Float64, maxElementStepSize::Float64)
        this = new()

        (elementIndex < 1) && throw(ArgumentError("Element index must be positive"))
        (sign(initialStepSize) == sign(maxElementStepSize)) || throw(ArgumentError("Step sizes must have same sign"))
        this.elementName = elementName
        this.elementIndex = Int16(elementIndex)
        this.initialStepSize = initialStepSize
        this.minStepSize = (initialStepSize < 0) ? -1E-10 : 1E-10
        this.maxStepSize = 1E-1
        this.maxElementStepSize = maxElementStepSize
        this.scaleFactor = 2.0
        this.minIterations = Int16(12)
        this.maxIterations = Int16(3)

        return this
    end
end
Base.:(==)(adaptiveStepSizeByElementGenerator1::AdaptiveStepSizeByElementGenerator, adaptiveStepSizeByElementGenerator2::AdaptiveStepSizeByElementGenerator) = ((adaptiveStepSizeByElementGenerator1.elementIndex == adaptiveStepSizeByElementGenerator2.elementIndex) && (adaptiveStepSizeByElementGenerator1.elementName == adaptiveStepSizeByElementGenerator2.elementName) && (adaptiveStepSizeByElementGenerator1.initialStepSize == adaptiveStepSizeByElementGenerator2.initialStepSize) && (adaptiveStepSizeByElementGenerator1.maxElementStepSize == adaptiveStepSizeByElementGenerator2.maxElementStepSize) && (adaptiveStepSizeByElementGenerator1.maxIterations == adaptiveStepSizeByElementGenerator2.maxIterations) && (adaptiveStepSizeByElementGenerator1.maxStepSize == adaptiveStepSizeByElementGenerator2.maxStepSize) && (adaptiveStepSizeByElementGenerator1.minIterations == adaptiveStepSizeByElementGenerator2.minIterations) && (adaptiveStepSizeByElementGenerator1.minStepSize == adaptiveStepSizeByElementGenerator2.minStepSize) && (adaptiveStepSizeByElementGenerator1.scaleFactor == adaptiveStepSizeByElementGenerator2.scaleFactor))

"""
    NumberStepsContinuationEndCheck(maxSteps)

Number of steps continuation end check object

# Arguments
- `maxSteps::Int64`: Maximum number of steps
"""
struct NumberStepsContinuationEndCheck <: AbstractContinuationEndCheck
    maxSteps::Int16                                                     # Maximum continuation steps

    function NumberStepsContinuationEndCheck(maxSteps::Int64)
        this = new(Int16(maxSteps))

        return this
    end
end
Base.:(==)(numberStepsContinuationEndCheck1::NumberStepsContinuationEndCheck, numberStepsContinuationEndCheck2::NumberStepsContinuationEndCheck) = (numberStepsContinuationEndCheck1.maxSteps == numberStepsContinuationEndCheck2.maxSteps)

"""
    BoundingBoxContinuationEndCheck(paramName, paramBounds)

Bounding box continuation end check object

# Arguments
- `paramName::String`: Parameter name
- `paramBounds::Matrix{Float64}`: Parameter lower/upper bounds
"""
mutable struct BoundingBoxContinuationEndCheck <: AbstractContinuationEndCheck
    paramBounds::Matrix{Float64}                                        # Parameter upper/lower bounds
    paramName::String                                                   # Parameter name
    variableBounds::Dict{Int16, Vector{Float64}}                        # Map between free variable index and bounds

    function BoundingBoxContinuationEndCheck(paramName::String, paramBounds::Matrix{Float64})
        this = new()

        this.paramName = paramName
        this.paramBounds = paramBounds
        this.variableBounds = Dict{Int16, Vector{Float64}}()

        return this
    end
end
Base.:(==)(boundingBoxContinuationEndCheck1::BoundingBoxContinuationEndCheck, boundingBoxContinuationEndCheck2::BoundingBoxContinuationEndCheck) = (isequal(boundingBoxContinuationEndCheck1.paramBounds, boundingBoxContinuationEndCheck2.paramBounds) && (boundingBoxContinuationEndCheck1.paramName == boundingBoxContinuationEndCheck2.paramName) && (boundingBoxContinuationEndCheck1.variableBounds == boundingBoxContinuationEndCheck2.variableBounds))

"""
    BundingBoxJumpCheck(paramName, paramBounds)

Bounding box jump check object

# Arguments
- `paramName::String`: Parameter name
- `paramBounds::Matrix{Float64}`: Parameter lower/upper bounds
"""
mutable struct BoundingBoxJumpCheck <: AbstractContinuationJumpCheck
    paramBounds::Matrix{Float64}                                        # Parameter upper/lower bounds
    paramName::String                                                   # Parameter name
    variableBounds::Dict{Int16, Vector{Float64}}                        # Map between free variable index and bounds

    function BoundingBoxJumpCheck(paramName::String, paramBounds::Matrix{Float64})
        this = new()

        this.paramName = paramName
        this.paramBounds = paramBounds
        this.variableBounds = Dict{Int16, Vector{Float64}}()

        return this
    end
end
Base.:(==)(boundingBoxJumpCheck1::BoundingBoxJumpCheck, boundingBoxJumpCheck2::BoundingBoxJumpCheck) = (isequal(boundingBoxJumpCheck1.paramBounds, boundingBoxJumpCheck2.paramBounds) && (boundingBoxJumpCheck1.paramName == boundingBoxJumpCheck2.paramName) && (boundingBoxJumpCheck1.variableBounds == boundingBoxJumpCheck2.variableBounds))

"""
    CR3BPNaturalParameterContinuationEngine(solution1, solution2, paramName, paramIndex, initialParamStepSize, maxParamStepSize; tol)

CR3BP natural parameter continuation engine object

# Arguments
- `solution1::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
- `solution2::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
- `paramName::String`: Natural parameter name
- `paramIndex::Int64`: Natural parameter index to step in
- `initialParamStepSize::Float64`: Initial step size
- `maxParamStepSize::Float64`: Maximum parameter step size
- `tol::Float64`: Convergence tolerance (default = 1E-10)
"""
mutable struct CR3BPNaturalParameterContinuationEngine
    corrector::CR3BPMultipleShooter                                     # Multiple shooter corrector for family
    dataInProgress::CR3BPContinuationData                               # Continuation data
    endChecks::Vector{AbstractContinuationEndCheck}                     # Continuation end checks
    jumpChecks::Vector{AbstractContinuationJumpCheck}                   # Continuation jump checks
    printProgress::Bool                                                 # Print progress?
    stepSizeGenerator::AdaptiveStepSizeByElementGenerator               # Step size generator
    storeIntermediateMembers::Bool                                      # Store intermediate family members?

    function CR3BPNaturalParameterContinuationEngine(solution1::CR3BPMultipleShooterProblem, solution2::CR3BPMultipleShooterProblem, paramName::String, paramIndex::Int64, initialParamStepSize::Float64, maxParamStepSize::Float64, tol::Float64 = 1E-10)
        this = new()

        this.corrector = CR3BPMultipleShooter(tol)
        this.dataInProgress = CR3BPContinuationData(solution1, solution2)
        this.stepSizeGenerator = AdaptiveStepSizeByElementGenerator(paramName, paramIndex, initialParamStepSize, maxParamStepSize)
        this.jumpChecks = []
        this.endChecks = []
        this.storeIntermediateMembers = true
        this.printProgress = true

        return this
    end
end
Base.:(==)(naturalParameterContinuationEngine1::CR3BPNaturalParameterContinuationEngine, naturalParameterContinuationEngine2::CR3BPNaturalParameterContinuationEngine) = ((naturalParameterContinuationEngine1.corrector == naturalParameterContinuationEngine2.corrector) && (naturalParameterContinuationEngine1.dataInProgress == naturalParameterContinuationEngine2.dataInProgress) && (naturalParameterContinuationEngine1.endChecks == naturalParameterContinuationEngine2.endChecks) && (naturalParameterContinuationEngine1.jumpChecks == naturalParameterContinuationEngine2.jumpChecks) && (naturalParameterContinuationEngine1.stepSizeGenerator == naturalParameterContinuationEngine2.stepSizeGenerator))

"""
    JacobiConstantContinuationEngine(solution1, solution2, initialParamStepSize, maxParamStepSize; tol)

Jacobi constant continuation engine object

# Arguments
- `solution1::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
- `solution2::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem solution
- `initialParamStepSize::Float64`: Initial step size
- `maxParamStepSize::Float64`: Maximum parameter step size
- `tol::Float64`: Convergence tolerance (default = 1E-10)
"""
mutable struct JacobiConstantContinuationEngine
    corrector::CR3BPMultipleShooter                                     # Multiple shooter corrector for family
    dataInProgress::CR3BPContinuationData                               # Continuation data
    endChecks::Vector{AbstractContinuationEndCheck}                     # Continuation end checks
    jumpChecks::Vector{AbstractContinuationJumpCheck}                   # Continuation jump checks
    printProgress::Bool                                                 # Print progress?
    stepSizeGenerator::AdaptiveStepSizeByElementGenerator               # Step size generator
    storeIntermediateMembers::Bool                                      # Store intermediate family members?

    function JacobiConstantContinuationEngine(solution1::CR3BPMultipleShooterProblem, solution2::CR3BPMultipleShooterProblem, initialParamStepSize::Float64, maxParamStepSize::Float64, tol::Float64 = 1E-10)
        this = new()

        this.corrector = CR3BPMultipleShooter(tol)
        this.dataInProgress = CR3BPContinuationData(solution1, solution2)
        this.stepSizeGenerator = AdaptiveStepSizeByElementGenerator("Jacobi Constant", 1, initialParamStepSize, maxParamStepSize)
        this.jumpChecks = []
        this.endChecks = []
        this.storeIntermediateMembers = true
        this.printProgress = true

        return this
    end
end
Base.:(==)(jacobiConstantContinuationEngine1::JacobiConstantContinuationEngine, jacobiConstantContinuationEngine2::JacobiConstantContinuationEngine) = ((jacobiConstantContinuationEngine1.corrector == jacobiConstantContinuationEngine2.corrector) && (jacobiConstantContinuationEngine1.dataInProgress == jacobiConstantContinuationEngine2.dataInProgress) && (jacobiConstantContinuationEngine1.endChecks == jacobiConstantContinuationEngine2.endChecks) && (jacobiConstantContinuationEngine1.jumpChecks == jacobiConstantContinuationEngine2.jumpChecks) && (jacobiConstantContinuationEngine1.stepSizeGenerator == jacobiConstantContinuationEngine2.stepSizeGenerator))

# """
#     Bifurcation(family, orbit, index, type, bifurcation)

# Bifurcation object

# # Arguments
# - `family::AbstractStructureFamily`: Structure family object
# - `orbit::AbstractTrajectoryStructure`: Trajectory structure object
# - `index::Int64`: Orbit index
# - `type::BifurcationType`: Bifurcation type
# - `bifurcation::Int64`: Bifurcation identifier
# """
# mutable struct Bifurcation
#     family::AbstractStructureFamily                         # Original family
#     FVVStep::Vector{Float64}                                # Free variable vector step into new family
#     ICStep::Vector{Float64}                                 # Initial condition step into new family
#     number::Int64                                           # Bifurcation identifier
#     orbit::AbstractTrajectoryStructure                      # Bifurcating structure
#     sortedEigenvalues::Vector{Complex{Float64}}             # Family-sorted eigenvalues
#     sortedEigenvectors::Matrix{Complex{Float64}}            # Family-sorted eigenvectors
#     type::BifurcationType                                   # Bifurcation type
    
#     function Bifurcation(family::AbstractStructureFamily, orbit::AbstractTrajectoryStructure, index::Int64, type::BifurcationType, bifurcation::Int64)
#         return new(family, Vector{Float64}(undef, getNumberFreeVariables(orbit.problem)), Vector{Float64}(undef, 6), bifurcation, orbit, family.eigenvalues[index], family.eigenvectors[index], type)
#     end
# end

"""
    CR3BPPeriodicOrbit(dynamicsModel, initialCondition, period, monodromy)

CR3BP periodic orbit object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `initialCondition::Vector{Float64}`: Initial condition [ndim]
- `period::Float64`: Period [ndim]
- `monodromy::Matrix{Float64}`: Monodromy matrix [ndim]
"""
struct CR3BPPeriodicOrbit
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    initialCondition::Vector{Float64}                                   # Initial condition [ndim]
    monodromy::StaticArrays.SMatrix{6, 6, Float64}                      # Monodromy matrix [ndim]
    period::Float64                                                     # Period [ndim]

    function CR3BPPeriodicOrbit(dynamicsModel::CR3BPDynamicsModel, initialCondition::Vector{Float64}, period::Float64, monodromy::Matrix{Float64})
        this = new(dynamicsModel, copy(initialCondition), StaticArrays.SMatrix{6, 6, Float64}(monodromy), copy(period))

        return this
    end
end
Base.:(==)(periodicOrbit1::CR3BPPeriodicOrbit, periodicOrbit2::CR3BPPeriodicOrbit) = ((periodicOrbit1.dynamicsModel == periodicOrbit2.dynamicsModel) && (periodicOrbit1.initialCondition == periodicOrbit2.initialCondition) && (periodicOrbit1.monodromy == periodicOrbit2.monodromy) && (periodicOrbit1.period == periodicOrbit2.period))

"""
    CR3BPOrbitFamily(dynamicsModel)

CR3BP orbit family object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
mutable struct CR3BPOrbitFamily
    # bifurcations::Vector{Bifurcation}                                   # Bifurcations
    dynamicsModel::CR3BPDynamicsModel                                   # CR3BP dynamics model object
    eigenvalues::Vector{Vector{Complex{Float64}}}                       # Eigenvalues [ndim]
    eigenvectors::Vector{Matrix{Complex{Float64}}}                      # Eigenvectors [ndim]
    hasBeenSorted::Bool                                                 # Eigendata has been sorted?
    initialConditions::Vector{Vector{Float64}}                          # Initial conditions [ndim]
    monodromies::Vector{Matrix{Float64}}                                # Monodromy matrices [ndim]
    periods::Vector{Float64}                                            # Periods [ndim]

    function CR3BPOrbitFamily(dynamicsModel::CR3BPDynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
        this.initialConditions = []
        this.periods = []
        this.monodromies = []
        this.eigenvalues = []
        this.eigenvectors = []
        this.hasBeenSorted = false
        # this.bifurcations = []

        return this
    end
end
Base.:(==)(orbitFamily1::CR3BPOrbitFamily, orbitFamily2::CR3BPOrbitFamily) = ((orbitFamily1.dynamicsModel == orbitFamily2.dynamicsModel) && (orbitFamily1.initialConditions == orbitFamily2.initialConditions)  && (orbitFamily1.monodromies == orbitFamily2.monodromies) && (orbitFamily1.periods == orbitFamily2.periods))

"""
    CR3BPManifoldArc(periodicOrbit, orbitTime, d, initialCondition; TOF)

CR3BP manifold arc object

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: Underlying CR3BP periodic orbit
- `orbitTime::Float64`: Time along orbit from initial condition [ndim]
- `d::Float64`: Step-off distance [ndim]
- `initialCondition::Vector{Complex{Float64}}`: Initial conditions [ndim]
- `TOF::Float64`: Time-of-flight [ndim] (default = 0.0)
"""
mutable struct CR3BPManifoldArc
    d::Float64                                                          # Step-off distance [ndim]
    initialCondition::Vector{Complex{Float64}}                          # Initial conditions [ndim]
    orbitTime::Float64                                                  # Normalized time along orbit from initial condition
    periodicOrbit::CR3BPPeriodicOrbit                                   # Underlying periodic orbit
    TOF::Float64                                                        # Time-of-flight [ndim]

    function CR3BPManifoldArc(periodicOrbit::CR3BPPeriodicOrbit, orbitTime::Float64, d::Float64, initialCondition::Vector{Complex{Float64}}, TOF::Float64 = 0.0)
        this = new()

        this.periodicOrbit = periodicOrbit
        this.orbitTime = orbitTime
        this.d = d
        this.initialCondition = initialCondition
        this.TOF = TOF

        return this
    end
end

"""
    CR3BPManifold(periodicOrbit, stability, direction; TOF)

CR3BP manifold object

# Arguments
- `periodicOrbit::CR3BPPeriodicOrbit`: Underlying CR3BP periodic orbit
- `stability::String`: Stability
- `direction::String`: Step-off direction
- `TOF::Float64`: Time-of-flight [ndim] (default = 0.0)
"""
mutable struct CR3BPManifold
    direction::String                                                   # Step-off direction
    ds::Vector{Float64}                                                 # Step-off distances [ndim]
    initialConditions::Vector{Vector{Complex{Float64}}}                 # Initial conditions [ndim]
    orbitTimes::Vector{Float64}                                         # Normalized times along orbit from initial condition
    periodicOrbit::CR3BPPeriodicOrbit                                   # Underlying periodic orbit
    stability::String                                                   # Stability
    TOF::Float64                                                        # Time-of-flight [ndim]

    function CR3BPManifold(periodicOrbit::CR3BPPeriodicOrbit, stability::String, direction::String, TOF::Float64 = 0.0)
        this = new()

        this.periodicOrbit = periodicOrbit
        this.orbitTimes = []
        this.ds = []
        this.stability = stability
        this.direction = direction
        this.initialConditions = []
        this.TOF = TOF

        return this
    end
end

"""
    BCR4BPSystemData(p1, p2, p4, b1)

BCR4BP system data object

# Arguments
- `p1::String`: Name of first primary
- `p2::String`: Name of second primary
- `p4::String`: Name of fourth primary
- `b1::String`: Name of first barycenter
"""
mutable struct BCR4BPSystemData
    primaryData::StaticArrays.SVector{4, BodyData}                      # Primary data objects
    primaryNames::StaticArrays.SVector{4, String}                       # Primary names
    primarySpiceIDs::StaticArrays.SVector{4, Int16}                     # Primary SPICE IDs

    function BCR4BPSystemData(p1::String, p2::String, p4::String, b1::String)
        this = new()

        this.primaryData = StaticArrays.SVector(BodyData(p1), BodyData(p2), BodyData(p4), BodyData(b1))
        this.primaryNames = StaticArrays.SVector(p1, p2, p4, b1)
        this.primarySpiceIDs = StaticArrays.SVector(this.primaryData[1].spiceID, this.primaryData[2].spiceID, this.primaryData[3].spiceID, this.primaryData[4].spiceID)
        (this.primaryData[2].parentSpiceID == this.primarySpiceIDs[1]) || throw(ArgumentError("First primary must be parent of second primary"))
        (this.primaryData[1].parentSpiceID == this.primarySpiceIDs[3]) || throw(ArgumentError("Fourth primary must be parent of fourth primary"))
        (this.primaryData[4].parentSpiceID == this.primarySpiceIDs[3]) || throw(ArgumentError("Fourth primary must be parent of first barycenter"))

        return this
    end
end
Base.:(==)(systemData1::BCR4BPSystemData, systemData2::BCR4BPSystemData) = ((systemData1.primaryData == systemData2.primaryData) && (systemData1.primaryNames == systemData2.primaryNames) && (systemData1.primarySpiceIDs == systemData2.primarySpiceIDs))

"""
    BCR4BP12DynamicsModel(systemData)

BCR4BP P1-P2 dynamics model object

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
struct BCR4BP12DynamicsModel
    systemData::BCR4BPSystemData                                        # BCR4BP system data object

    function BCR4BP12DynamicsModel(systemData::BCR4BPSystemData)
        this = new(systemData)

        return this
    end
end
Base.:(==)(dynamicsModel1::BCR4BP12DynamicsModel, dynamicsModel2::BCR4BP12DynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

"""
    BCR4BP12EquationsOfMotion(equationType, dynamicsModel)

BCR4BP P1-P2 EOM object

# Arguments
- `equationType::EquationType`: EOM type
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
struct BCR4BP12EquationsOfMotion
    dynamicsModel::BCR4BP12DynamicsModel                                # BCR4BP P1-P2 dynamics model object
    equationType::EquationType                                          # EOM type
    
    function BCR4BP12EquationsOfMotion(equationType::EquationType, dynamicsModel::BCR4BP12DynamicsModel)
        this = new(dynamicsModel, equationType)

        return this
    end
end
Base.:(==)(EOMs1::BCR4BP12EquationsOfMotion, EOMs2::BCR4BP12EquationsOfMotion) = ((EOMs1.dynamicsModel == EOMs2.dynamicsModel) && (EOMs1.equationType == EOMs2.equationType))

"""
    BCR4BP12Arc(dynamicsModel)

BCR4BP P1-P2 arc object

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
mutable struct BCR4BP12Arc
    dynamicsModel::BCR4BP12DynamicsModel                                # BCR4BP P1-P2 dynamics model object
    states::Vector{Vector{Float64}}                                     # State vectors along arc [ndim]
    times::Vector{Float64}                                              # Times along arc [ndim]

    function BCR4BP12Arc(dynamicsModel::BCR4BP12DynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
        this.states = []
        this.times = []

        return this
    end
end
Base.:(==)(arc1::BCR4BP12Arc, arc2::BCR4BP12Arc) = ((arc1.dynamicsModel == arc2.dynamicsModel) && (arc1.states == arc2.states) && (arc1.times == arc2.times))

"""
    BCR4BP41DynamicsModel(systemData)

BCR4BP P4-B1 dynamics model object

# Arguments
- `systemData::BCR4BPSystemData`: BCR4BP system data object
"""
struct BCR4BP41DynamicsModel
    systemData::BCR4BPSystemData                                        # BCR4BP system data object

    function BCR4BP41DynamicsModel(systemData::BCR4BPSystemData)
        this = new(systemData)

        return this
    end
end
Base.:(==)(dynamicsModel1::BCR4BP41DynamicsModel, dynamicsModel2::BCR4BP41DynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

"""
    BCR4BP41EquationsOfMotion(equationType, dynamicsModel)

BCR4BP P4-B1 EOM object

# Arguments
- `equationType::EquationType`: EOM type
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
struct BCR4BP41EquationsOfMotion
    dynamicsModel::BCR4BP41DynamicsModel                                # BCR4BP P4-B1 dynamics model object
    equationType::EquationType                                          # EOM type
    
    function BCR4BP41EquationsOfMotion(equationType::EquationType, dynamicsModel::BCR4BP41DynamicsModel)
        this = new(dynamicsModel, equationType)

        return this
    end
end
Base.:(==)(EOMs1::BCR4BP41EquationsOfMotion, EOMs2::BCR4BP41EquationsOfMotion) = ((EOMs1.dynamicsModel == EOMs2.dynamicsModel) && (EOMs1.equationType == EOMs2.equationType))

"""
    BCR4BP41Arc(dynamicsModel)

BCR4BP P4-B1 arc object

# Arguments
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
"""
mutable struct BCR4BP41Arc
    dynamicsModel::BCR4BP41DynamicsModel                                # BCR4BP P4-B1 dynamics model object
    states::Vector{Vector{Float64}}                                     # State vectors along arc [ndim]
    times::Vector{Float64}                                              # Times along arc [ndim]

    function BCR4BP41Arc(dynamicsModel::BCR4BP41DynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
        this.states = []
        this.times = []

        return this
    end
end
Base.:(==)(arc1::BCR4BP41Arc, arc2::BCR4BP41Arc) = ((arc1.dynamicsModel == arc2.dynamicsModel) && (arc1.states == arc2.states) && (arc1.times == arc2.times))

# """
#     TBPSystemData(p)

# TBP system object

# # Arguments
# - `p::String`: Name of primary
# """
# mutable struct TBPSystemData <: AbstractSystemData
#     charLength::Float64                                     # Characteristic length [km]
#     charTime::Float64                                       # Characteristic time [s]
#     gravParam::Float64                                      # Gravitational parameter [km^3/s^2]
#     numPrimaries::Int64                                     # Number of primaries that must exist in this system
#     primaryName::String                                     # Primary name
#     primarySpiceID::Int64                                   # Primary SPICE ID

#     function TBPSystemData(p::String)
#         this = new()

#         pData = BodyData(p)
#         this.numPrimaries = 1
#         this.primaryName = pData.name
#         this.primarySpiceID = pData.spiceID
#         this.gravParam = pData.gravParam
#         this.charLength = pData.orbitRadius
#         this.charTime = sqrt(this.charLength^3/this.gravParam)
        
#         return this
#     end
# end
# Base.:(==)(systemData1::TBPSystemData, systemData2::TBPSystemData) = (systemData1.primarySpiceID == systemData2.primarySpiceID)

# """
#     TBPDynamicsModel(systemData)

# TBP dynamics model object

# # Arguments
# - `systemData::TBPSystemData`: TBP system object
# """
# struct TBPDynamicsModel <: AbstractDynamicsModel
#     systemData::TBPSystemData                               # CR3BP system object

#     function TBPDynamicsModel(systemData::TBPSystemData)
#         return new(systemData)
#     end
# end
# Base.:(==)(dynamicsModel1::TBPDynamicsModel, dynamicsModel2::TBPDynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

# """
#     TBPEquationsOfMotion(equationType, dynamicsModel)

# TBP EOM object

# # Arguments
# - `equationType::EquationType`: EOM type
# - `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
# """
# struct TBPEquationsOfMotion <: AbstractEquationsOfMotion
#     dim::Int64                                              # State vector dimension
#     equationType::EquationType                              # EOM type
#     mu::Float64                                             # TBP gravitational parameter

#     function TBPEquationsOfMotion(equationType::EquationType, dynamicsModel::TBPDynamicsModel)
#         return new(getStateSize(dynamicsModel, equationType), equationType, dynamicsModel.systemData.gravParam)
#     end
# end
# Base.:(==)(EOMs1::TBPEquationsOfMotion, EOMs2::TBPEquationsOfMotion) = ((EOMs1.equationType == EOMs2.equationType) && (EOMs1.mu == EOMs2.mu))

# """
#     TBPTrajectory(initialCondition, dynamicsModel)

# TBP trajectory object

# # Arguments
# - `initialCondition::Vector{Float64}`: Initial conditions
# - `dynamicsModel::TBPDynamicsModel`: Dynamics model object
# """
# mutable struct TBPTrajectory <: AbstractTrajectoryStructure
#     a::Float64                                              # Semimajor axis [km]
#     dynamicsModel::TBPDynamicsModel                         # TBP Dynamics model object
#     E::Float64                                              # Energy [km^2/s^2]
#     e::Float64                                              # Eccentricity
#     h::Float64                                              # Angular momentum [km^2/s]
#     i::Float64                                              # Inclination [rad]
#     initialCondition::Vector{Float64}                       # Initial conditions
#     Omega::Float64                                          # Longitude of ascending node [rad]
#     omega::Float64                                          # Argument of periapsis [rad]
#     r_a::Float64                                            # Radius of apoapse [km]
#     r_p::Float64                                            # Radius of periapse [km]
#     theta::Float64                                          # True anomaly [rad]
#     TOF::Float64                                            # Time of flight

#     function TBPTrajectory(initialCondition::Vector{Float64}, dynamicsModel::TBPDynamicsModel)
#         return new(0.0, dynamicsModel, 0.0, 0.0, 0.0, 0.0, initialCondition, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     end
# end
# Base.:(==)(trajectory1::TBPTrajectory, trajectory2::TBPTrajectory) = ((trajectory1.a == trajectory2.a) && (trajectory1.dynamicsModel == trajectory2.dynamicsModel) && (trajectory1.E == trajectory2.E) && (trajectory1.e == trajectory2.e) && (trajectory1.h == trajectory2.h) && (trajectory1.i == trajectory2.i) && (trajectory1.initialCondition == trajectory2.initialCondition) && (trajectory1.Omega == trajectory2.Omega) && (trajectory1.omega == trajectory2.omega) && (trajectory1.theta == trajectory2.theta))

# include("bifurcation/Bifurcation.jl")
include("BCR4BP12/Arc12.jl")
include("BCR4BP12/DynamicsModel12.jl")
include("BCR4BP12/EquationsOfMotion12.jl")
include("BCR4BP12/SystemData.jl")
include("BCR4BP41/Arc41.jl")
include("BCR4BP41/DynamicsModel41.jl")
include("BCR4BP41/EquationsOfMotion41.jl")
include("continuation/AdaptiveStepSizeByElementGenerator.jl")
include("continuation/BoundingBoxContinuationEndCheck.jl")
include("continuation/BoundingBoxJumpCheck.jl")
include("continuation/NumberStepsContinuationEndCheck.jl")
include("corrections/ConstraintVectorL2NormConvergenceCheck.jl")
include("corrections/LeastSquaresUpdateGenerator.jl")
include("corrections/MinimumNormUpdateGenerator.jl")
include("corrections/StateMatchConstraint.jl")
include("corrections/Variable.jl")
include("CR3BP/Arc.jl")
include("CR3BP/ContinuationData.jl")
include("CR3BP/ContinuationFamily.jl")
include("CR3BP/ContinuityConstraint.jl")
include("CR3BP/DynamicsModel.jl")
include("CR3BP/EquationsOfMotion.jl")
include("CR3BP/JacobiConstantContinuationEngine.jl")
include("CR3BP/JacobiConstraint.jl")
include("CR3BP/Manifold.jl")
include("CR3BP/ManifoldArc.jl")
include("CR3BP/MultipleShooter.jl")
include("CR3BP/MultipleShooterProblem.jl")
include("CR3BP/NaturalParameterContinuationEngine.jl")
include("CR3BP/Node.jl")
include("CR3BP/OrbitFamily.jl")
include("CR3BP/PeriodicOrbit.jl")
include("CR3BP/Segment.jl")
include("CR3BP/StateConstraint.jl")
include("CR3BP/SystemData.jl")
include("propagation/EventFunctions.jl")
include("propagation/Propagator.jl")
include("spice/BodyName.jl")
include("spice/SpiceFunctions.jl")
# include("TBP/DynamicsModel.jl")
# include("TBP/EquationsOfMotion.jl")
# include("TBP/Trajectory.jl")
include("utilities/UtilityFunctions.jl")

end # module MBD
