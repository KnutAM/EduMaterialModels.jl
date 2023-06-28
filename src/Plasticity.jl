@kwdef struct J2Plasticity{E,T} <: AbstractMaterial
    # Elasticity
    e::E    = LinearIsotropicElasticity(;E=200e3, ν=0.3)
    Y0::T   = 200.0 # Initial yield limit
    Hiso::T = 10e3  # Isotropic hardening modulus
    κ∞::T   = 100.0 # Isotropic saturation stress
    Hkin::T = 30e3  # Kinematic hardening modulus 
    β∞::T   = 200.0 # Kinematic saturation stress
end

struct J2PlasticityState{T}  <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6}
    β::SymmetricTensor{2,3,T,6}
    κ::T
end
function MMB.initial_material_state(::J2Plasticity)
    ϵp = zero(SymmetricTensor{2,3})
    β = zero(SymmetricTensor{2,3})
    return J2PlasticityState(ϵp, β, 0.0)
end

function to_unknowns(::J2Plasticity, Δλ, ϵp, κ, β)
    ϵp_s = tosmandel(ϵp)
    β_s = tosmandel(β)
    return vcat(vcat(Δλ, ϵp_s), vcat(κ, β_s))
end

function from_unknowns(::J2Plasticity, x::SVector)
    Δλ = x[1]
    ϵp = frommandel(SymmetricTensor{2,3}, x; offset=1)
    κ = x[8]
    β = frommandel(SymmetricTensor{2,3}, x; offset=8)
    return Δλ, ϵp, κ, β
end

function process_solution(m::J2Plasticity, ϵ, old, Δt, x, drdx)
    _, ϵp, κ, β = from_unknowns(m, x)
    σ = elastic_stress(m.e, ϵ - ϵp)
    function σ_fun(x)
        _, ϵp_inside, _, _ = from_unknowns(m, x)
        return tosmandel(elastic_stress(m.e, ϵ - ϵp_inside))
    end
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x*dxdϵ
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x dxdϵ
    ∂σ∂ϵ = elastic_stiffness(m.e)
    ∂σ∂x = ForwardDiff.jacobian(σ_fun, x)
    rf_ϵs(ϵs) = residual(m, frommandel(SymmetricTensor{2,3}, ϵs), old, Δt, x)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵs, tosmandel(ϵ))
    dσdϵ = ∂σ∂ϵ - frommandel(SymmetricTensor{4,3}, ∂σ∂x * (drdx \ ∂r∂ϵ))
    return σ, dσdϵ, J2PlasticityState(ϵp, β, κ)
end

function MMB.material_response(m::J2Plasticity, ϵ, old, Δt=nothing, cache=get_cache(m), extras=NoExtraOutput(); options=Dict())
    ϵe_trial = ϵ - old.ϵp
    σ_trial = elastic_stress(m.e, ϵe_trial)
    Φ_trial = von_mises(σ_trial - old.β) - (m.Y0 + old.κ)
    if Φ_trial < 0 # Elastic response 
        return σ_trial, elastic_stiffness(m.e, ϵe_trial), old 
    else # Plastic 
        rf(x) = residual(m, ϵ, old, Δt, x)
        x0 = to_unknowns(m, 0.0, old.ϵp, old.κ, old.β)
        x, drdx, converged = newtonsolve(x0, rf)
        converged || throw(MMB.ConvergenceError("J2Plasticity did not converge locally"))
        σ, dσdϵ, new_state =  process_solution(m, ϵ, old, Δt, x, drdx)
        return σ, dσdϵ, new_state
    end
end

function residual(m::J2Plasticity, ϵ, old, Δt, x::SVector)
    Δλ, ϵp, κ, β = from_unknowns(m, x)
    σ = elastic_stress(m.e, ϵ-ϵp)
    ν, σ_vm = gradient(von_mises, σ-β, :all)
    RΦ = σ_vm - (m.Y0 + κ)
    Rϵp = (ϵp - old.ϵp) - Δλ*ν
    Rκ = (κ - old.κ) - Δλ*m.Hiso*(1 - κ/m.κ∞)
    Rβ = (β - old.β) - Δλ*((3/2)*m.Hkin)*(ν - (2/3)*β/m.β∞)
    return to_unknowns(m, RΦ, Rϵp, Rκ, Rβ)
end

