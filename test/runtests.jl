using EduMaterialModels
using Test, Tensors
import EduMaterialModels as EMM
import MaterialModelsBase as MMB

@testset "EduMaterialModels.jl" begin
    E = 200.e3
    m = EMM.J2Plasticity(;e=EMM.LinearIsotropicElasticity(;E=E, ν=0.3))
    ϵ = SymmetricTensor{2,1}.(tuple.(range(0, 0.002, 100)))
    stress_state = MMB.UniaxialStress()
    t = range(0, 1.0, length(ϵ))
    σ = EMM.simulate_response(m, stress_state, ϵ, t)
    @test first(σ[2]-σ[1]) / first(ϵ[2]-ϵ[1]) ≈ E

    # Test Explicit implementation
    mw = EMM.ExplicitPlasticity(m)
    σw = EMM.simulate_response(mw, stress_state, ϵ, t)
    @test all(isapprox.(σ, σw; atol=0.1)) # Explicit is bad ...

    # Test Rate-dependent implementation 
    m_rdep = EMM.J2ViscoPlasticity(;
        e=m.e, Y0=m.Y0, Hiso=m.Hiso, κ∞=m.κ∞, Hkin=m.Hkin, β∞=m.β∞, # as above
        tstar=1e-3, n=2.0)
    σ_fast = EMM.simulate_response(m_rdep, stress_state, ϵ, t)
    σ_slow = EMM.simulate_response(m_rdep, stress_state, ϵ, t*1e6)
    @test all(isapprox.(σ, σ_slow; atol=1e-3))
    @test all(s1[1] >= s2[1] for (s1, s2) in zip(σ_fast, σ_slow))
    @test all(s1[1] >= s2[1] for (s1, s2) in zip(σ_slow, σ))
end
