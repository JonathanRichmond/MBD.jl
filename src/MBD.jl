"""
Multi-Body Dynamics astrodynamics package

Author: Jonathan Richmond
C: 9/1/22
U: 7/29/23
"""
module MBD

import Base: ==
import DifferentialEquations, LightXML

const GRAVITY = 6.67384e-20
const UNINITIALIZED_INDEX = 0

"""
Enumerated type for EOMs
"""
@enum EquationType FULL SIMPLE STM

"""
Enumerated type for integrators
"""
@enum IntegratorType AB AM BS DP5 DP8

"""
Abstract type for constraints
"""
abstract type AbstractConstraint end

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
Abstract type for nonlinear problems
"""
abstract type AbstractNonlinearProblem end

"""
Abstract type for system objects
"""
abstract type AbstractSystemData end

"""
    BodyName(name)

Body name object

# Arguments
- `name::String`: Name of body
"""
struct BodyName
    name::String                                        # Name

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
    bodyRadius::Float64                                 # Body mean radius [km]
    gravParam::Float64                                  # Gravitational parameter [kg^3/s^2]
    inc::Float64                                        # Orbit inclination relative to Ecliptic J2000 frame [rad]
    mass::Float64                                       # Mass [kg]
    name::String                                        # Name
    orbitRadius::Float64                                # Orbit mean radius [km]
    parentSpiceID::Int64                                # Parent body SPICE ID
    RAAN::Float64                                       # Orbit Right-Ascension of Ascending Node relative to Ecliptic J2000 frame [rad]
    spiceID::Int64                                      # SPICE ID

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
    charLength::Float64                                 # Characteristic length [km]
    charMass::Float64                                   # Characteristic mass [kg]
    charTime::Float64                                   # Characteritic time [s]
    numPrimaries::Int64                                 # Number of primaries that must exist in this system
    params::Vector{Float64}                             # Other system parameters
    primaryNames::Vector{String}                        # Primary names
    primarySpiceIDs::Vector{Int64}                      # Primary SPICE IDs

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
    systemData::CR3BPSystemData                         # CR3BP system object

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
    dim::Int64                                          # State vector dimension
    equationType::EquationType                          # EOM type
    mu::Float64                                         # CR3BP system mass ratio

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
    integrator                                          # Integrator object
    integratorType::IntegratorType                      # Integrator type
    numSteps::Int64                                     # Number of steps for multi-step method

    function IntegratorFactory()
        return new(DifferentialEquations.DP8(), DP8, 3)
    end
end

"""
    Propagator()

Propagator object
"""
mutable struct Propagator
    absTol::Float64                                     # Absolute tolerance
    equationType::EquationType                          # EOM type
    events::Vector{AbstractEvent}                       # Integration events
    integratorFactory::IntegratorFactory                # Integrator object
    maxEvaluationCount::Int64                           # Maximum number of equation evaluations
    maxStep::Int64                                      # Maximum step size
    minEventTime::Float64                               # Minimum time between initial time and first event occurrence
    relTol::Float64                                     # Relative tolerance

    function Propagator()
        return new(1E-14, SIMPLE, [], IntegratorFactory(), typemax(Int64), 100, 1E-12, 1E-12)
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
    dynamicsModel::AbstractDynamicsModel                # Dynamics model object
    params::Vector{Float64}                             # System parameters
    states::Vector{Vector{Float64}}                     # State vectors along arc [ndim]
    times::Vector{Float64}                              # Times along arc [ndim]

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
    data::Vector{Float64}                               # Data values
    freeVarMask::Vector{Bool}                           # Free variable mask
    name::String                                        # Name

    function Variable(data::Vector{Float64}, freeVarMask::Vector{Bool})
        (length(freeVarMask) == length(data)) || throw(ArgumentError("Free variable mask length, $(length(freeVarMask)), must match data values length, $(length(data))"))

        return new(copy(data), copy(freeVarMask), "")
    end
end
Base.:(==)(variable1::Variable, variable2::Variable) = ((variable1.data == variable2.data) && (variable1.freeVarMask == variable2.freeVarMask))

"""
    Node(epoch, state, dynamicsModel)

Node object

# Arguments
- `epoch::Float64`: Epoch [ndim]
- `state::Vector{Float64}`: State vector [ndim]
- `dynamicsModel::AbstractDynamicsModel`: Dynamics model object
"""
mutable struct Node
    dynamicsModel::AbstractDynamicsModel                # Dynamics model object
    epoch::Variable                                     # Epoch
    state::Variable                                     # State

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
mutable struct Segment
    originNode::Node                                    # Origin node object
    propArc::Arc                                        # Propagated arc object
    propagator::Propagator                              # Propagator object
    propParams::Variable                                # Propagation parameters
    terminalNode::Node                                  # Terminal node object
    TOF::Variable                                       # Time of flight

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
    constraintIndexMap::Dict{AbstractConstraint, Int64} # Map between constraints and first index of contraint equation in constraint vector
    constraintVector::Vector{Float64}                   # Constraints
    freeVariableIndexMap::Dict{Variable, Int64}         # Map between free variables and first index of free variables in free variable vector
    freeVariableVector::Vector{Float64}                 # Free variables
    hasBeenBuilt::Bool                                  # Has been built?
    jacobian::Vector{Vector{Float64}}                   # Jacobian matrix
    nodes::Vector{Node}                                 # Nodes
    segments::Vector{Segment}                           # Segments

    function MultipleShooterProblem()
        return new(Dict{AbstractConstraint, Int64}(), [], Dict{Variable, Int64}(), [], false, [[]], [], [])
    end
end

include("corrections/MultipleShooterProblem.jl")
include("corrections/Node.jl")
include("corrections/Segment.jl")
include("corrections/Variable.jl")
include("CR3BP/DynamicsModel.jl")
include("CR3BP/EquationsOfMotion.jl")
include("CR3BP/SystemData.jl")
include("propagation/Arc.jl")
include("propagation/Propagator.jl")
include("spice/BodyName.jl")
include("utilities/UtilityFunctions.jl")

end # module MBD
