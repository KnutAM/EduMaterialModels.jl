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
    mw = EMM.ExplicitWrapper(m)
    σw = EMM.simulate_response(mw, stress_state, ϵ, t)
    @test all(isapprox.(σ, σw; atol=0.1)) # Explicit is bad ...
end
