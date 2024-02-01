"""
Multi-Body Dynamics astrodynamics package

Author: Jonathan Richmond
C: 9/1/22
U: 7/29/23
"""
module MBD

import Base: ==
import Combinatorics, DifferentialEquations, LightXML, LinearAlgebra, SPICE

const GRAVITY = 6.67384E-20
const UNINITIALIZED_INDEX = 0

"""
Enumerated type for bifurcations
"""
@enum BifurcationType TANGENT DOUBLING TRIPLING QUADRUPLING HOPF

"""
Enumerated type for EOMs
"""
@enum EquationType ARCLENGTH FULL SIMPLE STM

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
Abstract type for continuation engines
"""
abstract type AbstractContinuationEngine end

"""
Abstract type for continuation jump checks
"""
abstract type AbstractContinuationJumpCheck end

"""
Abstract type for convergence checks
"""
abstract type AbstractConvergenceCheck end

"""
Abstract type for dynamics models
"""
abstract type AbstractDynamicsModel end

"""
Abstract type for EOMs
"""
abstract type AbstractEquationsOfMotion end

"""
Abstract type for events
"""
abstract type AbstractEvent end

"""
Abstract type for linear solution generators
"""
abstract type AbstractLinearSolutionGenerator end

"""
Abstract type for nonlinear problems
"""
abstract type AbstractNonlinearProblem end

"""
Abstract type for nonlinear problem solvers
"""
abstract type AbstractNonlinearProblemSolver end

"""
Abstract type for structure families
"""
abstract type AbstractStructureFamily end

"""
Abstract type for system objects
"""
abstract type AbstractSystemData end

"""
Abstract type for targeters
"""
abstract type AbstractTargeter end

"""
Abstract type for trajectory structures
"""
abstract type AbstractTrajectoryStructure end

"""
Abstract type for update generators
"""
abstract type AbstractUpdateGenerator end

"""
Abstract type for objects containing variables
"""
abstract type IHasVariables end

"""
    BodyName(name)

Body name object

# Arguments
- `name::String`: Name of body
"""
struct BodyName
    name::String                                            # Name

    function BodyName(name)
        return new(name)
    end
end

"""
    BodyData(name)

Body object

# Arguments
- `name::String`: Name of body
"""
mutable struct BodyData
    bodyRadius::Float64                                     # Body mean radius [km]
    gravParam::Float64                                      # Gravitational parameter [kg^3/s^2]
    inc::Float64                                            # Orbit inclination relative to Ecliptic J2000 frame [rad]
    mass::Float64                                           # Mass [kg]
    name::String                                            # Name
    orbitRadius::Float64                                    # Orbit mean radius [km]
    parentSpiceID::Int64                                    # Parent body SPICE ID
    RAAN::Float64                                           # Orbit Right-Ascension of Ascending Node relative to Ecliptic J2000 frame [rad]
    spiceID::Int64                                          # SPICE ID

    function BodyData(name::String)
        this = new()

        dataFile::String = joinpath(@__DIR__,"body_data.xml")
        bodyName = BodyName(name)
        id::Int64 = getIDCode(bodyName)
        try
            doc::LightXML.XMLDocument = LightXML.parse_file(dataFile)
            root::LightXML.XMLElement = LightXML.root(doc)
            nodeList::Vector{LightXML.XMLElement} = LightXML.get_elements_by_tagname(root, "body")
            for node::LightXML.XMLElement in nodeList
                if LightXML.is_elementnode(node)
                    elID::String = LightXML.content(LightXML.get_elements_by_tagname(node, "id")[1])
                    this.spiceID = parse(Int64, elID)
                    if this.spiceID == id
                        this.name = LightXML.content(LightXML.get_elements_by_tagname(node, "name")[1])
                        parentIDStr::String = LightXML.content(LightXML.get_elements_by_tagname(node, "parentId")[1])
                        (parentIDStr == "NaN") || (this.parentSpiceID = parse(Int64, parentIDStr))
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

"""
    CR3BPSystemData(p1, p2)

CR3BP system object

# Arguments
- `p1::String`: Name of first primary
- `p2::String`: Name of second primary
"""
mutable struct CR3BPSystemData <: AbstractSystemData
    charLength::Float64                                     # Characteristic length [km]
    charMass::Float64                                       # Characteristic mass [kg]
    charTime::Float64                                       # Characteritic time [s]
    numPrimaries::Int64                                     # Number of primaries that must exist in this system
    params::Vector{Float64}                                 # Other system parameters
    primaryNames::Vector{String}                            # Primary names
    primarySpiceIDs::Vector{Int64}                          # Primary SPICE IDs

    function CR3BPSystemData(p1::String, p2::String)
        this = new()

        p1Data = BodyData(p1)
        p2Data = BodyData(p2)
        this.numPrimaries = 2
        this.primaryNames = [p1Data.name, p2Data.name]
        this.primarySpiceIDs = [p1Data.spiceID, p2Data.spiceID]
        (p2Data.parentSpiceID == this.primarySpiceIDs[1]) || throw(ArgumentError("First primary must be parent of second primary"))
        totalGravParam::Float64 = p1Data.gravParam+p2Data.gravParam
        this.charLength = p2Data.orbitRadius
        this.charMass = totalGravParam/GRAVITY
        this.charTime = sqrt(this.charLength^3/totalGravParam)
        this.params = [p2Data.gravParam/totalGravParam]
        
        return this
    end
end
Base.:(==)(systemData1::CR3BPSystemData, systemData2::CR3BPSystemData) = (systemData1.primarySpiceIDs == systemData2.primarySpiceIDs)

"""
    CR3BPDynamicsModel(systemData)

CR3BP dynamics model object

# Arguments
- `systemData::CR3BPSystemData`: CR3BP system object
"""
struct CR3BPDynamicsModel <: AbstractDynamicsModel
    systemData::CR3BPSystemData                             # CR3BP system object

    function CR3BPDynamicsModel(systemData::CR3BPSystemData)
        return new(systemData)
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
struct CR3BPEquationsOfMotion <: AbstractEquationsOfMotion
    dim::Int64                                              # State vector dimension
    equationType::EquationType                              # EOM type
    mu::Float64                                             # CR3BP system mass ratio

    function CR3BPEquationsOfMotion(equationType::EquationType, dynamicsModel::CR3BPDynamicsModel)
        return new(getStateSize(dynamicsModel, equationType), equationType, getMassRatio(dynamicsModel.systemData))
    end
end
Base.:(==)(EOMs1::CR3BPEquationsOfMotion, EOMs2::CR3BPEquationsOfMotion) = ((EOMs1.equationType == EOMs2.equationType) && (EOMs1.mu == EOMs2.mu))

"""
    IntegratorFactory()

Integrator object
"""
mutable struct IntegratorFactory
    integrator                                              # Integrator object
    integratorType::IntegratorType                          # Integrator type
    numSteps::Int64                                         # Number of steps for multi-step method

    function IntegratorFactory()
        return new(DifferentialEquations.DP8(), DP8, 3)
    end
end

"""
    Propagator()

Propagator object
"""
mutable struct Propagator
    absTol::Float64                                         # Absolute tolerance
    equationType::EquationType                              # EOM type
    events::Vector{AbstractEvent}                           # Integration events
    integratorFactory::IntegratorFactory                    # Integrator object
    maxEvaluationCount::Int64                               # Maximum number of equation evaluations
    maxStep::Int64                                          # Maximum step size
    minEventTime::Float64                                   # Minimum time between initial time and first event occurrence
    relTol::Float64                                         # Relative tolerance

    function Propagator()
        return new(1E-12, SIMPLE, [], IntegratorFactory(), typemax(Int64), 100, 1E-12, 1E-12)
    end
end
Base.:(==)(propagator1::Propagator, propagator2::Propagator) = ((propagator1.equationType == propagator2.equationType) && (propagator1.events == propagator2.events))

"""
    Arc(dynamicsModel)

Arc object

# Arguments
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
"""
mutable struct Arc
    dynamicsModel::AbstractDynamicsModel                    # Dynamics model object
    params::Vector{Float64}                                 # System parameters
    states::Vector{Vector{Float64}}                         # State vectors along arc [ndim]
    times::Vector{Float64}                                  # Times along arc [ndim]

    function Arc(dynamicsModel::AbstractDynamicsModel)
        return new(dynamicsModel, [], [[]], [])
    end
end
Base.:(==)(arc1::Arc, arc2::Arc) = ((arc1.dynamicsModel == arc2.dynamicsModel) && (arc1.params == arc2.params) && (arc1.states == arc2.states) && (arc1.times == arc2.times))

"""
    Variable(data, freeVarMask)

Variable object

# Arguments
- `data::Vector{Float64}`: Data values
- `freeVarMask::Vector{Bool}`: Free variable mask
"""
mutable struct Variable
    data::Vector{Float64}                                   # Data values
    freeVariableMask::Vector{Bool}                          # Free variable mask
    name::String                                            # Name

    function Variable(data::Vector{Float64}, freeVariableMask::Vector{Bool})
        (length(freeVariableMask) == length(data)) || throw(ArgumentError("Free variable mask length, $(length(freeVariableMask)), must match data values length, $(length(data))"))

        return new(copy(data), copy(freeVariableMask), "")
    end
end
Base.:(==)(variable1::Variable, variable2::Variable) = ((variable1.data == variable2.data) && (variable1.freeVariableMask == variable2.freeVariableMask))

"""
    Node(epoch, state, dynamicsModel)

Node object

# Arguments
- `epoch::Float64`: Epoch [ndim]
- `state::Vector{Float64}`: State vector [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
"""
mutable struct Node <: IHasVariables
    dynamicsModel::AbstractDynamicsModel                    # Dynamics model object
    epoch::Variable                                         # Epoch
    state::Variable                                         # State

    function Node(epoch::Float64, state::Vector{Float64}, dynamicsModel::AbstractDynamicsModel)
        this = new()

        this.dynamicsModel = dynamicsModel
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

        return this
    end
end
Base.:(==)(node1::Node, node2::Node) = ((node1.dynamicsModel == node2.dynamicsModel) && (node1.epoch == node2.epoch) && (node1.state == node2.state))

"""
    Segment(TOF, originNode, terminalNode)

Segment object

# Arguments
- `TOF::Float64`: Time-of-flight
- `originNode::Node`: Origin node
- `terminalNode::Node`: Terminal node
"""
mutable struct Segment <: IHasVariables
    originNode::Node                                        # Origin node object
    propArc::Arc                                            # Propagated arc object
    propagator::Propagator                                  # Propagator object
    propParams::Variable                                    # Propagation parameters
    terminalNode::Node                                      # Terminal node object
    TOF::Variable                                           # Time of flight

    function Segment(TOF::Float64, originNode::Node, terminalNode::Node)
        this = new()

        this.TOF = Variable([TOF], [true])
        this.TOF.name = "Segment Time-of-Flight"
        this.propParams = Variable(Array{Float64}(undef, 0), Array{Bool}(undef, 0))
        this.propParams.name = "Segment Propagation Parameters"
        this.propArc = Arc(originNode.dynamicsModel)
        (originNode == terminalNode) && throw(ArgumentError("Origin and terminal nodes cannot be identical"))
        this.originNode = originNode
        this.terminalNode = terminalNode
        this.propagator = Propagator()
        this.propagator.equationType = SIMPLE

        return this
    end
end
Base.:(==)(segment1::Segment, segment2::Segment) = ((segment1.originNode == segment2.originNode) && (segment1.terminalNode == segment2.terminalNode) && (segment1.TOF == segment2.TOF))

"""
    MultipleShooterProblem()

Multiple shooter problem object
"""
mutable struct MultipleShooterProblem <: AbstractNonlinearProblem
    constraintIndexMap::Dict{AbstractConstraint, Int64}     # Map between constraints and first index of contraint equation in constraint vector
    constraintVector::Vector{Float64}                       # Constraints
    freeVariableIndexMap::Dict{Variable, Int64}             # Map between free variables and first index of free variables in free variable vector
    freeVariableVector::Vector{Float64}                     # Free variables
    hasBeenBuilt::Bool                                      # Has been built?
    jacobian::Vector{Vector{Float64}}                       # Jacobian matrix
    nodes::Vector{Node}                                     # Nodes
    segments::Vector{Segment}                               # Segments

    function MultipleShooterProblem()
        return new(Dict{AbstractConstraint, Int64}(), [], Dict{Variable, Int64}(), [], false, [[]], [], [])
    end
end

"""
    ContinuityConstraint(segment)

Continuity constraint object

# Arguments
- `segment::Segment`: Segment object
"""
mutable struct ContinuityConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int64}                       # Constrained state indices
    segment::Segment                                        # Segment object

    function ContinuityConstraint(segment::Segment)
        return new(1:length(segment.originNode.state.data), segment)
    end
end
Base.:(==)(continuityConstraint1::ContinuityConstraint, continuityConstraint2::ContinuityConstraint) = ((continuityConstraint1.constrainedIndices == continuityConstraint2.constrainedIndices) && (continuityConstraint1.segment == continuityConstraint2.segment))

"""
    StateConstraint(node, indices, values)

State constraint object

# Arguments
- `node::Node`: Node object
- `indices::Vector{Int64}`: Constrained state indices
- `values::Vector{Float64}`: Constraint values
"""
mutable struct StateConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int64}                       # Constrained state indices
    values::Vector{Float64}                                 # Constraint values
    variable::Variable                                      # Constrained variable

    function StateConstraint(node::Node, indices::Vector{Int64}, values::Vector{Float64})
        this = new()

        this.variable = node.state
        (length(indices) == length(values)) || throw(ArgumentError("Number of indices, $(length(indices)), must match number of values, $(length(values))"))
        checkIndices(indices, length(this.variable.data))
        this.constrainedIndices = copy(indices)
        this.values = copy(values)
        freeVariableMask::Vector{Bool} = getFreeVariableMask(this.variable)
        for index::Int64 in this.constrainedIndices
            freeVariableMask[index] || throw(ArgumentError("Variable element $index is constrained but not free to vary"))
        end

        return this
    end
end
Base.:(==)(stateConstraint1::StateConstraint, stateConstraint2::StateConstraint) = ((stateConstraint1.constrainedIndices == stateConstraint2.constrainedIndices) && (stateConstraint1.values == stateConstraint2.values) && (stateConstraint1.variable == stateConstraint2.variable))

"""
    StateMatchConstraint(state1, state2, indices)

State match constraint object

# Arguments
- `state1::Variable`: First state vector
- `state2::Variable`: Second state vector
- `indices::Vector{Int64}`: Constrained state indices
"""
mutable struct StateMatchConstraint <: AbstractConstraint
    constrainedIndices::Vector{Int64}                       # Constrained state indices
    variable1::Variable                                     # First state variable
    variable2::Variable                                     # Second state variable

    function StateMatchConstraint(state1::Variable, state2::Variable, indices::Vector{Int64})
        this = new(indices, state1, state2)

        (length(this.variable1.data) == length(this.variable2.data)) || throw(ArgumentError("First state length, $(length(this.variable1.data)), must match second state length, $(length(this.variable2.data))"))
        checkIndices(this.constrainedIndices, length(this.variable1.data))
        mask1::Vector{Bool} = getFreeVariableMask(this.variable1)
        mask2::Vector{Bool} = getFreeVariableMask(this.variable2)
        for i::Int64 in eachindex(this.constrainedIndices)
            (mask1[this.constrainedIndices[i]] || mask2[this.constrainedIndices[i]]) || throw(ArgumentError("Cannot constrain index = $(this.constrainedIndices[i]); it is not a free variable"))
        end

        return this
    end
end
Base.:(==)(stateMatchConstraint1::StateMatchConstraint, stateMatchConstraint2::StateMatchConstraint) = ((stateMatchConstraint1.variable1 == stateMatchConstraint2.variable1) && (stateMatchConstraint1.variable2 == stateMatchConstraint2.variable2) && (stateMatchConstraint1.constrainedIndices == stateMatchConstraint2.constrainedIndices))

"""
    JacobiConstraint(node, value)

Jacobi constraint object

# Arguments
- `node::Node`: Node object
- `value::Float64`: Constraint value
"""
mutable struct JacobiConstraint <: AbstractConstraint
    dynamicsModel::CR3BPDynamicsModel                       # Dynamics model object
    epoch::Variable                                         # Epoch
    state::Variable                                         # Constrained state
    value::Float64                                          # Constraint value

    function JacobiConstraint(node::Node, value::Float64)
        this = new()

        this.epoch = node.epoch
        this.state = node.state
        (typeof(node.dynamicsModel) == CR3BPDynamicsModel) ? (this.dynamicsModel = node.dynamicsModel) : throw(ArgumentError("Expecting CR3BP dynamics model object"))
        this.value = value

        return this
    end
end
Base.:(==)(jacobiConstraint1::JacobiConstraint, jacobiConstraint2::JacobiConstraint) = ((jacobiConstraint1.dynamicsModel == jacobiConstraint2.dynamicsModel) && (jacobiConstraint1.state == jacobiConstraint2.state) && (jacobiConstraint1.epoch == jacobiConstraint2.epoch) && (jacobiConstraint1.value == jacobiConstraint2.value))

"""
    ConstraintVectorL2NormConvergenceCheck(tol)

Constraint vector L2 norm convergence check object

# Arguments
- `tol::Float64`: Convergence tolerance (optional)
"""
mutable struct ConstraintVectorL2NormConvergenceCheck <: AbstractConvergenceCheck
    maxVectorNorm::Float64                                  # Maximum allowable vector norm

    function ConstraintVectorL2NormConvergenceCheck(tol::Float64 = 1E-10)
        return new(tol)
    end
end

"""
    MinimumNormUpdateGenerator()

Minimum norm update generator object
"""
struct MinimumNormUpdateGenerator <: AbstractUpdateGenerator
    function MinimumNormUpdateGenerator()
        return new()
    end
end

"""
    LeastSquaresUpdateGenerator()

Least squares update generator object
"""
struct LeastSquaresUpdateGenerator <: AbstractUpdateGenerator
    function LeastSquaresUpdateGenerator()
        return new()
    end
end

"""
    MultipleShooter(tol)

Multiple shooter object

# Arguments
- `tol::Float64`: Convergence tolerance (optional)
"""
mutable struct MultipleShooter <: AbstractNonlinearProblemSolver
    convergenceCheck::AbstractConvergenceCheck              # Convergence check object
    maxIterations::Int64                                    # Maximum number of solver iterations
    printProgress::Bool                                     # Print progress?
    recentIterationCount::Int64                             # Number of iterations during last solve
    solutionInProgress::MultipleShooterProblem              # Multiple shooter problem being solved
    updateGenerators::Vector{AbstractUpdateGenerator}       # Update generators

    function MultipleShooter(tol::Float64 = 1E-10)
        return new(ConstraintVectorL2NormConvergenceCheck(tol), 25, false, 0, MultipleShooterProblem(), [MinimumNormUpdateGenerator(), LeastSquaresUpdateGenerator()])
    end
end

"""
    ContinuationData()

Continuation data object
"""
mutable struct ContinuationData
    converging::Bool                                        # Converging?
    currentStepSize::Float64                                # Current continuation step size
    familyMembers::Vector{MultipleShooterProblem}           # Computed family family members
    forceEndContinuation::Bool                              # Force end of continuation?
    fullStep::Vector{Float64}                               # Full step along free variable vector
    initialGuess::MultipleShooterProblem                    # First member passed to corrections algorithm
    nextGuess::MultipleShooterProblem                       # Inital guess for next family member
    numIterations::Int64                                    # Number of iterations required to converge previous solution
    previousSolution::MultipleShooterProblem                # Most recently converged family member
    stepCount::Int64                                        # Number of steps from initial guess
    twoPreviousSolution::MultipleShooterProblem             # Second previously converged family member

    function ContinuationData()
        return new(true, 1, [], false, [], MultipleShooterProblem(), MultipleShooterProblem(), 0, MultipleShooterProblem(), 0, MultipleShooterProblem())
    end
end

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
    elementIndex::Int64                                     # Free variable index
    elementName::String                                     # Free variable name
    initialStepSize::Float64                                # Initial continuation step size
    maxElementStepSize::Float64                             # Maximum continuation step size for element
    maxIterations::Int64                                    # Maximum iterations to converge solution
    maxStepSize::Float64                                    # Maximum continuation step size
    minIterations::Int64                                    # Minimum iterations to converge solution
    minStepSize::Float64                                    # Minimum continuation step size
    scaleFactor::Float64                                    # Step size scale scaleFactor

    function AdaptiveStepSizeByElementGenerator(elementName::String, elementIndex::Int64, initialStepSize::Float64, maxElementStepSize::Float64)
        this = new()

        this.initialStepSize = initialStepSize
        this.minStepSize = (initialStepSize < 0) ? -1E-10 : 1E-10
        this.maxStepSize = 1E-1
        this.scaleFactor = 2
        this.minIterations = 10
        this.maxIterations = 3
        this.elementName = elementName
        (elementIndex < 1) && throw(ArgumentError("Element index must be positive"))
        this.elementIndex = elementIndex
        this.maxElementStepSize = maxElementStepSize

        return this
    end
end

"""
    NumberStepsContinuationEndCheck(maxSteps)

Number of steps continuation end check object

# Arguments
- `maxSteps::Int64`: Maximum number of steps
"""
struct NumberStepsContinuationEndCheck <: AbstractContinuationEndCheck
    maxSteps::Int64                                         # Maximum continuation steps

    function NumberStepsContinuationEndCheck(maxSteps::Int64)
        return new(maxSteps)
    end
end

"""
    BoundingBoxContinuationEndCheck(paramName, paramBounds)

Bounding box continuation end check object

# Arguments
- `paramName::String`: Parameter name
- `paramBounds::Matrix{Float64}`: Parameter lower/upper bounds
"""
mutable struct BoundingBoxContinuationEndCheck <: AbstractContinuationEndCheck
    paramBounds::Matrix{Float64}                            # Parameter upper/lower bounds
    paramName::String                                       # Parameter name
    variableBounds::Dict{Int64, Vector{Float64}}            # Map between free variable index and bounds

    function BoundingBoxContinuationEndCheck(paramName::String, paramBounds::Matrix{Float64})
        return new(paramBounds, paramName, Dict{Int64, Vector{Float64}}())
    end
end

"""
    BundingBoxJumpCheck(paramName, paramBounds)

Bounding box jump check object

# Arguments
- `paramName::String`: Parameter name
- `paramBounds::Matrix{Float64}`: Parameter lower/upper bounds
"""
mutable struct BoundingBoxJumpCheck <: AbstractContinuationJumpCheck
    paramBounds::Matrix{Float64}                            # Parameter upper/lower bounds
    paramName::String                                       # Parameter name
    variableBounds::Dict{Int64, Vector{Float64}}            # Map between free variable index and bounds

    function BoundingBoxJumpCheck(paramName::String, paramBounds::Matrix{Float64})
        return new(paramBounds, paramName, Dict{Int64, Vector{Float64}}())
    end
end

"""
    NaturalParameterContinuationEngine(paramName, paramIndex, initialParamStepSize, maxParamStepSize; tol)

Natural parameter continuation engine object

# Arguments
- `paramName::String`: Natural parameter name
- `paramIndex::Int64`: Natural parameter index to step in
- `initialParamStepSize::Float64`: Initial step size
- `maxParamStepSize::Float64`: Maximum parameter step size
- `tol::Float64`: Convergence tolerance (optional)
"""
mutable struct NaturalParameterContinuationEngine <: AbstractContinuationEngine
    corrector::MultipleShooter                              # Multiple shooter corrector for family
    dataInProgress::ContinuationData                        # Continuation data
    endChecks::Vector{AbstractContinuationEndCheck}         # Continuation end checks
    jumpChecks::Vector{AbstractContinuationJumpCheck}       # Continuation jump checks
    printProgress::Bool                                     # Print progress?
    stepSizeGenerator::AdaptiveStepSizeByElementGenerator   # Step size generator
    storeIntermediateMembers::Bool                          # Store intermediate family members?

    function NaturalParameterContinuationEngine(paramName::String, paramIndex::Int64, initialParamStepSize::Float64, maxParamStepSize::Float64, tol::Float64 = 1E-10)
        return new(MultipleShooter(tol), ContinuationData(), [], [], true, AdaptiveStepSizeByElementGenerator(paramName, paramIndex, initialParamStepSize, maxParamStepSize), true)
    end
end

"""
    JacobiConstantContinuationEngine(initialParamStepSize, maxParamStepSize; tol)

Jacobi constant continuation engine object

# Arguments
- `initialParamStepSize::Float64`: Initial step size
- `maxParamStepSize::Float64`: Maximum parameter step size
- `tol::Float64`: Convergence tolerance (optional)
"""
mutable struct JacobiConstantContinuationEngine <: AbstractContinuationEngine
    corrector::MultipleShooter                              # Multiple shooter corrector for family
    dataInProgress::ContinuationData                        # Continuation data
    endChecks::Vector{AbstractContinuationEndCheck}         # Continuation end checks
    jumpChecks::Vector{AbstractContinuationJumpCheck}       # Continuation jump checks
    printProgress::Bool                                     # Print progress?
    stepSizeGenerator::AdaptiveStepSizeByElementGenerator   # Step size generator
    storeIntermediateMembers::Bool                          # Store intermediate family members?

    function JacobiConstantContinuationEngine(initialParamStepSize::Float64, maxParamStepSize::Float64, tol::Float64 = 1E-10)
        return new(MultipleShooter(tol), ContinuationData(), [], [], true, AdaptiveStepSizeByElementGenerator("Jacobi Constant", 1, initialParamStepSize, maxParamStepSize), true)
    end
end

"""
    CR3BPPeriodicOrbit(multipleShooterProblem, targeter)

CR3BP periodic orbit object

# Arguments
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
- `targeter::AbstractTargeter`: Targeter object
"""
mutable struct CR3BPPeriodicOrbit <: AbstractTrajectoryStructure
    BrouckeStability::Vector{Float64}                       # Broucke stability parameters
    eigenvalues::Vector{Complex{Float64}}                   # Monodromy matrix eigenvalues [ndim]
    eigenvectors::Matrix{Complex{Float64}}                  # Monodromy matrix eigenvectors [ndim]
    initialCondition::Vector{Float64}                       # Initial conditions [ndim]
    JacobiConstant::Float64                                 # Jacobi constant
    monodromy::Matrix{Float64}                              # Monodromy matrix [ndim]
    nu::Float64                                             # Maximum stability index [ndim]
    period::Float64                                         # Period [ndim]
    problem::MultipleShooterProblem                         # Solved multiple shooter problem
    targeter::AbstractTargeter                              # Targeter
    tau::Float64                                            # Time constant [ndim]

    function CR3BPPeriodicOrbit(multipleShooterProblem::MultipleShooterProblem, targeter::AbstractTargeter)
        return new(zeros(Float64, 2), Vector{Complex{Float64}}(undef, 6), Matrix{Complex{Float64}}(undef, 6, 6), Vector{Float64}(undef, 6), 0.0, Matrix{Float64}(undef, 6, 6), 0.0, 0.0, multipleShooterProblem, targeter, 0.0)
    end
end

"""
    CR3BPOrbitFamily(familyMembers)

CR3BP orbit family object

# Arguments
- `familyMembers::Vector{CR3BPPeriodicOrbit}`: Family members
"""
mutable struct CR3BPOrbitFamily <: AbstractStructureFamily
    alternateIndices::Vector{Vector{Complex{Float64}}}      # Alternate stability indices
    bifurcations::Vector{AbstractBifurcation}               # Bifurcations
    eigenvalues::Vector{Vector{Complex{Float64}}}           # Sorted family eigenvalues
    eigenvectors::Vector{Matrix{Complex{Float64}}}          # Sorted family eigenvectors
    familyMembers::Vector{CR3BPPeriodicOrbit}               # Family members
    stabilityIndices::Vector{Vector{Float64}}               # Stability indices

    function CR3BPOrbitFamily(familyMembers::Vector{CR3BPPeriodicOrbit})
        return new(Vector{Vector{Complex{Float64}}}(undef, length(familyMembers)), [], Vector{Vector{Complex{Float64}}}(undef, length(familyMembers)), Vector{Matrix{Complex{Float64}}}(undef, length(familyMembers)), familyMembers, Vector{Vector{Float64}}(undef, length(familyMembers)))
    end
end

"""
    CR3BPManifoldArc(initialCondition, periodicOrbit)

CR3BP manifold arc object

# Arguments
- `initialCondition::Vector{Complex{Float64}}`: Initial conditions [ndim]
- `periodicOrbit::CR3BPPeriodicOrbit`: Underlying CR3BP periodic orbit
"""
mutable struct CR3BPManifoldArc <: AbstractTrajectoryStructure
    initialCondition::Vector{Float64}                       # Initial conditions [ndim]
    JacobiConstant::Float64                                 # Jacobi constant
    periodicOrbit::CR3BPPeriodicOrbit                       # Underlying periodic orbit
    TOF::Float64                                            # Time of flight

    function CR3BPManifoldArc(initialCondition::Vector{Complex{Float64}}, periodicOrbit::CR3BPPeriodicOrbit)
        return new(real.(initialCondition), getJacobiConstant(periodicOrbit.targeter.dynamicsModel, real.(initialCondition)), periodicOrbit, 0.0)
    end
end

"""
    TBPSystemData(p)

TBP system object

# Arguments
- `p::String`: Name of primary
"""
mutable struct TBPSystemData <: AbstractSystemData
    charLength::Float64                                     # Characteristic length [km]
    charTime::Float64                                       # Characteristic time [s]
    gravParam::Float64                                      # Gravitational parameter [km^3/s^2]
    numPrimaries::Int64                                     # Number of primaries that must exist in this system
    primaryName::String                                     # Primary name
    primarySpiceID::Int64                                   # Primary SPICE ID

    function TBPSystemData(p::String)
        this = new()

        pData = BodyData(p)
        this.numPrimaries = 1
        this.primaryName = pData.name
        this.primarySpiceID = pData.spiceID
        this.gravParam = pData.gravParam
        this.charLength = pData.orbitRadius
        this.charTime = sqrt(this.charLength^3/this.gravParam)
        
        return this
    end
end
Base.:(==)(systemData1::TBPSystemData, systemData2::TBPSystemData) = (systemData1.primarySpiceID == systemData2.primarySpiceID)

"""
    TBPDynamicsModel(systemData)

TBP dynamics model object

# Arguments
- `systemData::TBPSystemData`: TBP system object
"""
struct TBPDynamicsModel <: AbstractDynamicsModel
    systemData::TBPSystemData                               # CR3BP system object

    function TBPDynamicsModel(systemData::TBPSystemData)
        return new(systemData)
    end
end
Base.:(==)(dynamicsModel1::TBPDynamicsModel, dynamicsModel2::TBPDynamicsModel) = (dynamicsModel1.systemData == dynamicsModel2.systemData)

"""
    TBPEquationsOfMotion(equationType, dynamicsModel)

TBP EOM object

# Arguments
- `equationType::EquationType`: EOM type
- `dynamicsModel::TBPDynamicsModel`: TBP dynamics model object
"""
struct TBPEquationsOfMotion <: AbstractEquationsOfMotion
    dim::Int64                                              # State vector dimension
    equationType::EquationType                              # EOM type
    mu::Float64                                             # TBP gravitational parameter

    function TBPEquationsOfMotion(equationType::EquationType, dynamicsModel::TBPDynamicsModel)
        return new(getStateSize(dynamicsModel, equationType), equationType, dynamicsModel.systemData.gravParam)
    end
end
Base.:(==)(EOMs1::TBPEquationsOfMotion, EOMs2::TBPEquationsOfMotion) = ((EOMs1.equationType == EOMs2.equationType) && (EOMs1.mu == EOMs2.mu))

"""
    TBPTrajectory(initialCondition, dynamicsModel)

TBP trajectory object

# Arguments
- `initialCondition::Vector{Float64}`: Initial conditions
- `dynamicsModel::TBPDynamicsModel`: Dynamics model object
"""
mutable struct TBPTrajectory <: AbstractTrajectoryStructure
    a::Float64                                              # Semimajor axis [km]
    dynamicsModel::TBPDynamicsModel                         # TBP Dynamics model object
    E::Float64                                              # Energy [km^2/s^2]
    e::Float64                                              # Eccentricity
    h::Float64                                              # Angular momentum [km^2/s]
    i::Float64                                              # Inclination [rad]
    initialCondition::Vector{Float64}                       # Initial conditions
    Omega::Float64                                          # Longitude of ascending node [rad]
    omega::Float64                                          # Argument of periapsis [rad]
    r_a::Float64                                            # Radius of apoapse [km]
    r_p::Float64                                            # Radius of periapse [km]
    theta::Float64                                          # True anomaly [rad]
    TOF::Float64                                            # Time of flight

    function TBPTrajectory(initialCondition::Vector{Float64}, dynamicsModel::TBPDynamicsModel)
        return new(0.0, dynamicsModel, 0.0, 0.0, 0.0, 0.0, initialCondition, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

"""
    Bifurcation(family, orbit, index, type, bifurcation)

Bifurcation object

# Arguments
- `family::AbstractStructureFamily`: Structure family object
- `orbit::AbstractTrajectoryStructure`: Trajectory structure object
- `index::Int64`: Orbit index
- `type::BifurcationType`: Bifurcation type
- `bifurcation::Int64`: Bifurcation identifier
"""
mutable struct Bifurcation
    family::AbstractStructureFamily                         # Original family
    FVVStep::Vector{Float64}                                # Free variable vector step into new family
    ICStep::Vector{Float64}                                 # Initial condition step into new family
    number::Int64                                           # Bifurcation identifier
    orbit::AbstractTrajectoryStructure                      # Bifurcating structure
    sortedEigenvalues::Vector{Complex{Float64}}             # Family-sorted eigenvalues
    sortedEigenvectors::Matrix{Complex{Float64}}            # Family-sorted eigenvectors
    type::BifurcationType                                   # Bifurcation type
    
    function Bifurcation(family::AbstractStructureFamily, orbit::AbstractTrajectoryStructure, index::Int64, type::BifurcationType, bifurcation::Int64)
        return new(family, Vector{Float64}(undef, getNumberFreeVariables(orbit.problem)), Vector{Float64}(undef, 6), bifurcation, orbit, family.eigenvalues[index], family.eigenvectors[index], type)
    end
end

include("bifurcation/Bifurcation.jl")
include("continuation/AdaptiveStepSizeByElementGenerator.jl")
include("continuation/BoundingBoxContinuationEndCheck.jl")
include("continuation/BoundingBoxJumpCheck.jl")
include("continuation/JacobiConstantContinuationEngine.jl")
include("continuation/NaturalParameterContinuationEngine.jl")
include("continuation/NumberStepsContinuationEndCheck.jl")
include("corrections/ConstraintVectorL2NormConvergenceCheck.jl")
include("corrections/ContinuityConstraint.jl")
include("corrections/LeastSquaresUpdateGenerator.jl")
include("corrections/MinimumNormUpdateGenerator.jl")
include("corrections/MultipleShooter.jl")
include("corrections/MultipleShooterProblem.jl")
include("corrections/Node.jl")
include("corrections/Segment.jl")
include("corrections/StateConstraint.jl")
include("corrections/StateMatchConstraint.jl")
include("corrections/Variable.jl")
include("CR3BP/DynamicsModel.jl")
include("CR3BP/EquationsOfMotion.jl")
include("CR3BP/JacobiConstraint.jl")
include("CR3BP/OrbitFamily.jl")
include("CR3BP/PeriodicOrbit.jl")
include("CR3BP/SystemData.jl")
include("propagation/Arc.jl")
include("propagation/EventFunctions.jl")
include("propagation/Propagator.jl")
include("spice/BodyName.jl")
include("spice/SpiceFunctions.jl")
include("TBP/DynamicsModel.jl")
include("TBP/EquationsOfMotion.jl")
include("TBP/Trajectory.jl")
include("utilities/UtilityFunctions.jl")

end # module MBD
