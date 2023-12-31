function von_mises(σ)
    σdev = dev(σ)
    return sqrt((3/2)*(σdev⊡σdev))
end

macaulay(x) = max(zero(x), x)