"""
EquationsOfMotion methods

Author: Jonathan LeFevre Richmond
C: 1/15/26
"""


function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{AbstractEquationsOfMotion, EquationType, Vararg}, t::Float64)
    Logging.@debug "computeDerivatives!(Abstract) called" type=typeof(params[1])

    Logging.@error "computeDerivatives! not implemented for abstract EquationsOfMotion" type=typeof(params[1])
    throw(ErrorException("computeDerivatives! not implemented for abstract EquationsOfMotion"))
end


function computeDerivatives!(qdot::Vector{Float64}, q::Vector{Float64}, params::Tuple{CR3BPEquationsOfMotion, EquationType, Vararg}, t::Float64)
    Logging.@debug "computeDerivatives!(CR3BP) called"

    if !(isa(qdot, AbstractVector) && all(isa.(qdot, Float64)))
        Logging.@error "qdot must be a Vector{Float64}" type=typeof(qdot)
        throw(ArgumentError("qdot must be a Vector{Float64}"))
    end
    if !(isa(q, AbstractVector) && all(isa.(q, Float64)))
        Logging.@error "q must be a Vector{Float64}" type=typeof(q)
        throw(ArgumentError("q must be a Vector{Float64}"))
    end
    if isempty(q)
        Logging.@error "q is empty"
        throw(ArgumentError("q must be non-empty"))
    end
    if isempty(qdot)
        Logging.@error "qdot is empty"
        throw(ArgumentError("qdot must be non-empty"))
    end
    if !all(isfinite, q)
        Logging.@error "State vector q contains non-finite values" has_nan=any(isnan, q) has_inf=any(isinf, q)
        throw(ArgumentError("q must contain only finite values"))
    end
    if !all(isfinite, qdot)
        Logging.@error "qdot contains non-finite values" has_nan=any(isnan, qdot) has_inf=any(isinf, qdot)
        throw(ArgumentError("qdot must contain only finite values"))
    end
    
    model::AbstractDynamicsModel = params[1].dynamicsModel
    equationType::EquationType = params[2]
    if model === nothing || !isa(model, CR3BPEquationsOfMotion)
        Logging.@error "Invalid CR3BPEquationsOfMotion supplied to computeDerivatives!" type=typeof(params[1])
        throw(ArgumentError("params[1] must be a valid CR3BPEquationsOfMotion instance"))
    end
    n_in_q::Int16 = Int16(length(q))
    n_in_qdot::Int16 = Int16(length(qdot))
    n_type::Int16 = Int16(getStateSize(model, equationType))
    if n_in_qdot != n_in_q
        Logging.@error "qdot and q vectors have mismatched sizes" qdot_length=length(qdot) q_length=length(q)
        throw(ArgumentError("qdot and q must have the same length"))
    end
    if n_in_qdot != n_type
        Logging.@error "Input state vector size does not match equation type size" n_in=n_in_qdot expected=n_type
        throw(ArgumentError("Input state vector size ($(n_in_qdot)) does not match equation type size ($(n_type))"))
    end

    try
        μ::Float64 = getMassRatios(model)
        if !isfinite(μ) || !(0 < μ < 1)
            Logging.@error "Invalid mass ratio computed" μ=μ
            throw(ErrorException("Mass ratio must be finite and in (0, 1)"))
        end
        r_13::Float64 = getDistance2Primary(model, 1, q[1:3])
        r_23::Float64 = getDistance2Primary(model, 2, q[1:3])
        if !isfinite(r_13) || r_13 <= 0
            Logging.@error "Invalid distance to primary 1" r_13=r_13
            throw(ErrorException("Distance to primary 1 must be positive and finite"))
        end
        if !isfinite(r_23) || r_23 <= 0
            Logging.@error "Invalid distance to primary 2" r_23=r_23
            throw(ErrorException("Distance to primary 2 must be positive and finite"))
        end
        r3_13::Float64 = r_13^3
        r3_23::Float64 = r_23^3

        # Compute SIMPLE state derivatives (always needed)
        qdot[1:3] = copy(q[4:6])
        qdot[4] = 2*q[5]+q[1]-(1-μ)*(q[1]+μ)/r3_13-μ*(q[1]+μ-1)/r3_23
        qdot[5] = q[2]-2*q[4]-(1-μ)*q[2]/r3_13-μ*q[2]/r3_23
        qdot[6] = -(1-μ)*q[3]/r3_13-μ*q[3]/r3_23
        if !all(isfinite, qdot[1:6])
            Logging.@error "Computed SIMPLE derivatives contain non-finite values" has_nan=any(isnan, qdot[1:6]) has_inf=any(isinf, qdot[1:6])
            throw(ErrorException("Computed SIMPLE derivatives must be finite"))
        end

        # Compute STM derivatives if needed
        if equationType != SIMPLE     
            r5_13::Float64 = r_13^5
            r5_23::Float64 = r_23^5
            pseudopotentialHessian::Matrix{Float64} = getPseudopotentialHessian(model, q[1:3])
            if !all(isfinite, pseudopotentialHessian)
                Logging.@error "Pseudopotential Hessian contains non-finite values"
                throw(ErrorException("Pseudopotential Hessian must contain only finite values"))
            end
            [qdot[6+6*(c-1)+r] = q[9+6*(c-1)+r] for r in 1:3, c in 1:6]
            for c::Int16 in 1:6
                qdot[10+6*(c-1)] = pseudopotentialHessian[1]*q[7+6*(c-1)]+pseudopotentialHessian[4]*q[8+6*(c-1)]+pseudopotentialHessian[5]*q[9+6*(c-1)]+2*q[11+6*(c-1)]
                qdot[11+6*(c-1)] = pseudopotentialHessian[4]*q[7+6*(c-1)]+pseudopotentialHessian[2]*q[8+6*(c-1)]+pseudopotentialHessian[6]*q[9+6*(c-1)]-2*q[10+6*(c-1)]
                qdot[12+6*(c-1)] = pseudopotentialHessian[5]*q[7+6*(c-1)]+pseudopotentialHessian[6]*q[8+6*(c-1)]+pseudopotentialHessian[3]*q[9+6*(c-1)]
            end
            if !all(isfinite, qdot[7:42])
                Logging.@error "Computed STM derivatives contain non-finite values" has_nan=any(isnan, qdot[7:42]) has_inf=any(isinf, qdot[7:42])
                throw(ErrorException("Computed STM derivatives must be finite"))
            end
        end

        # Compute additional derivatives based on equation type
        if equationType == ARCLENGTH
            qdot[43] = sqrt(q[4]^2+q[5]^2+q[6]^2)
            if !isfinite(qdot[43])
                Logging.@error "Arc-length derivative is non-finite" qdot_43=qdot[43]
                throw(ErrorException("Arc-length derivative must be finite"))
            end
        elseif equationType == MOMENTUM
            qdot[43] = q[1]*q[4]+q[2]*q[5]+q[3]*q[6]
            if !isfinite(qdot[43])
                Logging.@error "Momentum derivative is non-finite" qdot_43=qdot[43]
                throw(ErrorException("Momentum derivative must be finite"))
            end
        end
        Logging.@info "Successfully computed derivatives for CR3BP dynamics" equationType=equationType computed_elements=length(qdot)
    catch e
        Logging.@error "Failed to compute derivatives for CR3BP" exception=e
        rethrow()
    end
end
