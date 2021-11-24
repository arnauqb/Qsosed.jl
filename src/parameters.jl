export Parameters

struct Parameters{T<:AbstractFloat}
    # T float, S flag, U bool, V int, W string
    M::T
    mdot::T
    spin::T
end
function Parameters(;M, mdot, spin)
    return Parameters(M, mdot, spin)
end
