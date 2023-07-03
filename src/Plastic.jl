"""
    AbstractPlastic 

An `AbstractPlastic` materials have the same state variables,
`PlasticState`, unknowns for the local equation system `[Δλ, ϵp, κ, β]`,
and von mises yield criterion. They may have different evolution laws.
"""
abstract type AbstractPlastic  <: AbstractMaterial end

struct PlasticState{T}  <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6}
    β::SymmetricTensor{2,3,T,6}
    κ::T
end
function MMB.initial_material_state(::AbstractPlastic)
    ϵp = zero(SymmetricTensor{2,3})
    β = zero(SymmetricTensor{2,3})
    return PlasticState(ϵp, β, 0.0)
end

function to_unknowns(::AbstractPlastic, Δλ, ϵp, κ, β)
    ϵp_s = tosmandel(ϵp)
    β_s = tosmandel(β)
    return vcat(vcat(Δλ, ϵp_s), vcat(κ, β_s))
end

function from_unknowns(::AbstractPlastic, x::SVector)
    Δλ = x[1]
    ϵp = frommandel(SymmetricTensor{2,3}, x; offset=1)
    κ = x[8]
    β = frommandel(SymmetricTensor{2,3}, x; offset=8)
    return Δλ, ϵp, κ, β
end

function process_solution(m::AbstractPlastic, ϵ, old::PlasticState, Δt, x, drdx)
    _, ϵp, κ, β = from_unknowns(m, x)
    σ = elastic_stress(m, ϵ - ϵp)
    function σ_fun(x)
        _, ϵp_inside, _, _ = from_unknowns(m, x)
        return tosmandel(elastic_stress(m, ϵ - ϵp_inside))
    end
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x*dxdϵ
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x dxdϵ
    ∂σ∂ϵ = elastic_stiffness(m)
    ∂σ∂x = ForwardDiff.jacobian(σ_fun, x)
    rf_ϵs(ϵs) = residual(m, frommandel(SymmetricTensor{2,3}, ϵs), old, Δt, x)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵs, tosmandel(ϵ))
    dσdϵ = ∂σ∂ϵ - frommandel(SymmetricTensor{4,3}, ∂σ∂x * (drdx \ ∂r∂ϵ))
    return σ, dσdϵ, PlasticState(ϵp, β, κ)
end

function MMB.material_response(m::AbstractPlastic, ϵ, old::PlasticState, Δt=nothing, cache=get_cache(m), extras=NoExtraOutput(); options=Dict())
    ϵe_trial = ϵ - old.ϵp
    σ_trial = elastic_stress(m, ϵe_trial)
    Y0 = initial_yield_limit(m)
    Φ_trial = von_mises(σ_trial - old.β) - (Y0 + old.κ)
    if Φ_trial < 0 # Elastic response 
        return σ_trial, elastic_stiffness(m, ϵe_trial), old 
    else # Plastic 
        rf(x) = residual(m, ϵ, old, Δt, x)
        x0 = to_unknowns(m, 0.0, old.ϵp, old.κ, old.β)
        x, drdx, converged = newtonsolve(x0, rf)
        converged || throw(MMB.NoLocalConvergence("$(nameof(typeof(m))) did not converge locally"))
        σ, dσdϵ, new_state =  process_solution(m, ϵ, old, Δt, x, drdx)
        return σ, dσdϵ, new_state
    end
end

"""
    J2Plasticity(;e, Y0, Hiso, κ∞, Hkin, β∞)

The `J2Plasticity` material model assumes linear elasticity combined 
with rate independent von mises plasticity. Nonlinear isotropic (Voce)
and kinematic (Armstrong-Frederick) hardening  are included as well.
"""
@kwdef struct J2Plasticity{E,T} <: AbstractPlastic
    # Elasticity
    e::E    = LinearIsotropicElasticity(;E=200e3, ν=0.3)
    Y0::T   = 200.0 # Initial yield limit
    Hiso::T = 10e3  # Isotropic hardening modulus
    κ∞::T   = 100.0 # Isotropic saturation stress
    Hkin::T = 30e3  # Kinematic hardening modulus 
    β∞::T   = 200.0 # Kinematic saturation stress
end

elastic_stiffness(m::J2Plasticity, args...) = elastic_stiffness(m.e, args...)
elastic_stress(m::J2Plasticity, args...) = elastic_stress(m.e, args...)
initial_yield_limit(m::J2Plasticity) = m.Y0

function residual(m::J2Plasticity, ϵ, old::PlasticState, Δt, x::SVector)
    Δλ, ϵp, κ, β = from_unknowns(m, x)
    σ = elastic_stress(m.e, ϵ-ϵp)
    ν, σ_vm = gradient(von_mises, σ-β, :all)
    RΦ = σ_vm - (m.Y0 + κ)
    Rϵp = (ϵp - old.ϵp) - Δλ*ν
    Rκ = (κ - old.κ) - Δλ*m.Hiso*(1 - κ/m.κ∞)
    Rβ = (β - old.β) - Δλ*((2/3)*m.Hkin)*(ν - (3/2)*β/m.β∞)
    return to_unknowns(m, RΦ, Rϵp, Rκ, Rβ)
end

# Rate dependent response 
"""
    J2ViscoPlasticity(;e, Y0, Hiso, κ∞, Hkin, β∞, n, tstar)

The `J2ViscoPlasticity` material model assumes linear elasticity combined 
with rate dependent (Norton-type) von mises plasticity. Nonlinear isotropic (Voce)
and kinematic (Armstrong-Frederick) hardening  are included as well. Specifically,
the overstress function is 
```math
\\eta(\\Phi) = \\left[ \\frac{\\langle \\Phi \\rangle}{Y_0}\\right]^n
```
"""
@kwdef struct J2ViscoPlasticity{E,T} <: AbstractPlastic
    # Elasticity
    e::E    = LinearIsotropicElasticity(;E=200e3, ν=0.3)
    Y0::T    = 200.0 # Initial yield limit
    Hiso::T  = 10e3  # Isotropic hardening modulus
    κ∞::T    = 100.0 # Isotropic saturation stress
    Hkin::T  = 30e3  # Kinematic hardening modulus 
    β∞::T    = 200.0 # Kinematic saturation stress
    n::T     = 2.0   # Overstress exponent parameter
    tstar::T = 0.1   # Characteristic relaxation time
end

elastic_stiffness(m::J2ViscoPlasticity, args...) = elastic_stiffness(m.e, args...)
elastic_stress(m::J2ViscoPlasticity, args...) = elastic_stress(m.e, args...)
initial_yield_limit(m::J2ViscoPlasticity) = m.Y0

function residual(m::J2ViscoPlasticity, ϵ, old::PlasticState, Δt, x::SVector)
    Δλ, ϵp, κ, β = from_unknowns(m, x)
    σ = elastic_stress(m.e, ϵ-ϵp)
    ν, σ_vm = gradient(von_mises, σ-β, :all)
    Φ = σ_vm - (m.Y0 + κ)
    RΦ = Δλ - (Δt/m.tstar) * (macaulay(Φ)/m.Y0)^m.n
    Rϵp = (ϵp - old.ϵp) - Δλ*ν
    Rκ = (κ - old.κ) - Δλ*m.Hiso*(1 - κ/m.κ∞)
    Rβ = (β - old.β) - Δλ*((2/3)*m.Hkin)*(ν - (3/2)*β/m.β∞)
    return to_unknowns(m, RΦ, Rϵp, Rκ, Rβ)
end

