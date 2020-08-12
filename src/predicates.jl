function Base.hash(α::Cyclotomic, h::UInt)
    normalform!(α)
    return hash(coeffs(α), hash(conductor(α), hash(Cyclotomic, h)))
end

function Base.:(==)(α::Cyclotomic, β::Cyclotomic)
    coeffs(α) == coeffs(β) && return true

    if conductor(α) == conductor(β)
        normalform!(α)
        normalform!(β)
        return coeffs(α) == coeffs(β)
    else
        l = lcm(conductor(α), conductor(β))
        return embed(α, l) == embed(β, l)
    end
end

Base.:(==)(α::Cyclotomic{T,V}, x::R) where {T,V,R<:Real} =
    α == Cyclotomic(1, [x])
Base.:(==)(x::T, α::Cyclotomic) where {T<:Number} = α == x

function Base.isapprox(
    α::Cyclotomic{T},
    x::S;
    atol::Real = 0,
    rtol::Real = atol > 0 ? 0 : sqrt(max(eps(x), maximum(eps, coeffs(α)))),
) where {T<:AbstractFloat,S<:AbstractFloat}
    return isapprox(S(α), x; atol = atol, rtol = rtol)
end

function Base.isapprox(
    α::Cyclotomic{T},
    x::S;
    atol::Real = 0,
    rtol::Real = atol > 0 ? 0 : eps(x),
) where {T,S<:AbstractFloat}
    return isapprox(S(α), x; atol = atol, rtol = rtol)
end

Base.iszero(α::Cyclotomic) =
    all(iszero, values(α)) || (normalform!(α); all(iszero, values(α)))

Base.isreal(α::Cyclotomic) =
    α == conj(α) || conductor(reduced_embedding(α)) == 1

function Base.isone(α::Cyclotomic)
    β = reduced_embedding(α)
    conductor(β) == 1 || return false
    return isone(β[0])
end

"""
    isnormalized(α::Cyclotomic, basis)
Check if `α` is already in normal form with respect to the given basis.
"""
function isnormalized(α::Cyclotomic, basis = zumbroich(conductor(α)))
    # return all(in(basis), exponents(α))
    for e in exponents(α)
        e in basis || return false
    end
    return true
end

function _all_equal(α::Cyclotomic, exps, value = α[first(exps)])
    # return all(e -> α[e] == val, exps)
    for e in exps
        α[e] == value || return false
    end
    return true
end
