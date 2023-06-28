struct ExplicitWrapper{M<:AbstractMaterial} <: AbstractMaterial
    m::M
end
MMB.initial_material_state(m::ExplicitWrapper) = MMB.initial_material_state(m.m)

function MMB.material_response(mw::ExplicitWrapper{<:J2Plasticity}, ϵ, old, Δt=nothing, cache=get_cache(mw), extras=NoExtraOutput(); options=Dict())
    dσdϵ = gradient(e -> explicit_response(Val(false), mw.m, e, old), ϵ)
    σ, new = explicit_response(Val(true), mw.m, ϵ, old)
    return σ, dσdϵ, new
end

function explicit_response(return_state::Val, m::J2Plasticity, ϵ, old)
    E4 = elastic_stiffness(m.e)
    σ_trial = E4 ⊡ (ϵ - old.ϵp)
    Φ_trial = von_mises(σ_trial - old.β) - (m.Y0 + old.κ)
    if Φ_trial < 0 # Elastic response 
        isa(return_state, Val{false}) && return σ_trial
        return σ_trial, old
    else # Plastic
        ν = gradient(von_mises, σ_trial - old.β)

        # Evolution from the yield surface, goverened by dot(Φ) = 0
        dσdλ = - E4 ⊡ ν # dot(ϵp) = ν
        dκdλ = m.Hiso*(1 - old.κ/m.κ∞)
        dβdλ = (2/3)*m.Hkin*(ν - (3/2)*old.β/m.β∞)
        
        # dΦdλ = dΦdσ dσdλ + dΦdβ dβdλ + dΦdκ dκdλ # Full derivative for all states
        dΦdλ = ν ⊡ (dσdλ - dβdλ) - dκdλ
        # Φ_trial + dΦdλ * Δλ = 0
        Δλ = - Φ_trial/dΦdλ

        σ = σ_trial + Δλ * dσdλ
        isa(return_state, Val{false}) && return σ
        κ = old.κ + Δλ * dκdλ
        β = old.β + Δλ * dβdλ
        ϵp = old.ϵp + Δλ * ν
        return σ, J2PlasticityState(ϵp, β, κ)
    end
end