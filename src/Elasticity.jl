struct LinearIsotropicElasticity{T}
    G::T 
    K::T 
end

function LinearIsotropicElasticity(;kwargs...)
    length(kwargs) == 2 || throw(ArgumentError("Exactly two parameters must be given to LinearIsotropicElasticity"))
    G, K = transform_to_GK(kwargs, ;kwargs...)
    return LinearIsotropicElasticity(G, K)
end

function transform_to_GK(kwargs, ;G=nothing, K=nothing, E=nothing, ν=nothing)
    is_given(args...) = false
    is_given(::Number, ::Number) = true
    is_given(G, K) && return G, K 
    is_given(G, E) && return G, E*G/((3*(3*G-E)))
    is_given(E, ν) && return E/(2*(1+ν)), E/(3*(1-2*ν))
    is_given(K, E) && return 3K*E/(9K-E), K
    params = join([string(p) for (p,v) in kwargs], " and ")
    throw(ArgumentError("The input combination " * params * " is not supported"))
end

elastic_stress(m::LinearIsotropicElasticity, ϵ) = 2*m.G*dev(ϵ) + 3*m.K*vol(ϵ)
elastic_stiffness(m::LinearIsotropicElasticity, args...) = gradient(x->elastic_stress(m, x), zero(SymmetricTensor{2,3}))

