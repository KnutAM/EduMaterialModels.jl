# Static voigt format, see Tensors.jl #194
"""
    tosvoigt(A::Union{SecondOrderTensor,FourthOrderTensor}; offdiagscale=1)

`SVector` version of [`tovoigt`](@ref): converts the tensor `A` to a voigt `SVector` with the default voigt ordering. 
The keyword `offdiagscale` is only available for symmetric tensors. See also [`tosmandel`](@ref).
The output can be converted back to a tensor with the standard [`fromvoigt`](@ref) function.
"""
tosvoigt(A::Union{Tensor{2},Tensor{4}}) = _tosvoigt(A)
tosvoigt(A::Union{SymmetricTensor{2},SymmetricTensor{4}}; offdiagscale=one(eltype(A))) = _tosvoigt(A, offdiagscale)

"""
    tosmandel(A::Union{SecondOrderTensor,FourthOrderTensor})

`SVector` version of [`tomandel`](@ref): converts the tensor `A` to a mandel (`offdiagscale=√2` for symmetric tensors)
`SVector` with the default voigt ordering. See also [`tosvoigt`](@ref). 
The output can be converted back to a tensor with the standard [`frommandel`](@ref) function.
"""
tosmandel(A::Union{Tensor{2},Tensor{4}}) = _tosvoigt(A)
tosmandel(A::Union{SymmetricTensor{2},SymmetricTensor{4}}) = _tosvoigt(A, sqrt(2*one(eltype(A))))

@generated function _tosvoigt(A::TT, s::T=one(T)) where {TT<:SecondOrderTensor{dim,T}} where {dim,T}
    # Define internals for generation
    idx_fun(i, j) = Tensors.compute_index(Tensors.get_base(A), i, j)
    maxind(j) = TT<:SymmetricTensor ? j : dim
    N = Tensors.n_components(Tensors.get_base(A))

    exps = Expr(:tuple)
    append!(exps.args, [nothing for _ in 1:N]) # "Preallocate" to allow indexing directly

    for j in 1:dim, i in 1:maxind(j)
        voigt_ind = Tensors.DEFAULT_VOIGT_ORDER[dim][i,j]
        if i==j
            exps.args[voigt_ind] = :(Tensors.get_data(A)[$(idx_fun(i, j))])
        else
            exps.args[voigt_ind] = :(s*Tensors.get_data(A)[$(idx_fun(i, j))])
        end
    end

    quote
        $(Expr(:meta, :inline))
        @inbounds return SVector{$N, $T}($exps)
    end
end

@generated function _tosvoigt(A::TT, s::T=one(T)) where {TT<:FourthOrderTensor{dim,T}} where {dim,T}
    # Define internals for generation
    idx_fun(i, j, k, l) = Tensors.compute_index(Tensors.get_base(A), i, j, k, l)
    maxind(j) = TT<:SymmetricTensor ? j : dim
    voigt_lin_index(vi, vj) = (vj-1)*N + vi 
    N = Int(sqrt(Tensors.n_components(Tensors.get_base(A))))
    
    exps = Expr(:tuple)
    append!(exps.args, [nothing for _ in 1:N^2]) # "Preallocate" to allow indexing directly

    for l in 1:dim, k in 1:maxind(l), j in 1:dim, i in 1:maxind(j)
        voigt_lin_ind = voigt_lin_index(Tensors.DEFAULT_VOIGT_ORDER[dim][i,j], Tensors.DEFAULT_VOIGT_ORDER[dim][k,l])
        if i==j && k==l
            exps.args[voigt_lin_ind] = :(Tensors.get_data(A)[$(idx_fun(i, j, k, l))])
        elseif i!=j && k!=l
            exps.args[voigt_lin_ind] = :(s*s*Tensors.get_data(A)[$(idx_fun(i, j, k, l))])
        else
            exps.args[voigt_lin_ind] = :(s*Tensors.get_data(A)[$(idx_fun(i, j, k, l))])
        end
    end
    
    quote
        $(Expr(:meta, :inline))
        @inbounds return SMatrix{$N, $N, $T}($exps)
    end
end