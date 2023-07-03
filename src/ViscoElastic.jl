"""
    Zener1D(;E1, E2, η)

Special implementation of the Zener material model, which only works 
for uniaxial stress, i.e. when `stress_state=MaterialModelsBase.UniaxialStress()`
is given to `material_response`. 
"""
@kwdef struct Zener1D{T} <: AbstractMaterial
    E1::T = 20.e3
    E2::T = 10.e3
    η::T  = 0.1e3
end

struct Zener1DState{T}
    ϵv::T 
end
MaterialModelsBase.initial_material_state(::Zener1D{T}) where T = Zener1DState(zero(T))

viscous_strain(m::Zener1D, ϵ::Number, old::Zener1DState, Δt) = (Δt*m.E2*ϵ + m.η*old.ϵv)/(m.E2*Δt + m.η)

function calculate_stress(m::Zener1D, ϵt::SymmetricTensor{2,1}, old::Zener1DState, Δt)
    ϵ = first(ϵt)
    ϵv = viscous_strain(m, ϵ, old, Δt)
    σ = m.E1*ϵ + m.E2*(ϵ - ϵv)
    return SymmetricTensor{2,1}((σ,))
end

function MaterialModelsBase.material_response(::UniaxialStress, m::Zener1D, ϵ, old::Zener1DState, Δt, args...; kwargs...)
    dσdϵ, σ = gradient(ϵt -> calculate_stress(m, ϵt, old, Δt), ϵ, :all)
    new_state = Zener1DState(viscous_strain(m, first(ϵ), old, Δt))
    return σ, dσdϵ, new_state
end
