####
#   Arithmetic

zero!(α::Cyclotomic{T}) where {T} = (coeffs(α) .= zero(T); α)
function one!(α::Cyclotomic{T}) where {T}
    zero!(α)
    α[0] = one(α[0])
    return α
end
Base.zero(α::Cyclotomic) = zero!(similar(α))
Base.zero(α::Cyclotomic, m::Integer) = zero!(similar(α, m))
Base.one(α::Cyclotomic) = one!(similar(α))

############################
# Module structure:

Base.:-(α::Cyclotomic) = Cyclotomic(-coeffs(α))

for op in (:+, :-)
    @eval begin
        function Base.$op(α::Cyclotomic{T}, r::R) where {T,R<:Real}
            res = similar(α, promote_type(T, R))
            copyto!(coeffs(res), coeffs(α))
            res[0] = $op(res[0], r)
            return res
        end
    end
end

function Base.:-(r::R, α::Cyclotomic{T}) where {T,R<:Real}
    res = similar(α, promote_type(T, R))
    coeffs(res) .= -1 .* coeffs(α)
    res[0] += r
    return res
end

Base.:+(r::Real, α::Cyclotomic) = α + r

function mul!(out::Cyclotomic, α::Cyclotomic, c::Real)
    return (coeffs(out) .= coeffs(α) .* c; out)
end
function div!(out::Cyclotomic, α::Cyclotomic, c::Real)
    return (coeffs(out) .= div.(coeffs(α), c); out)
end

function Base.:*(c::T, α::Cyclotomic{S}) where {S,T<:Real}
    return mul!(similar(α, promote_type(S, T)), α, c)
end
Base.:*(α::Cyclotomic, c::T) where {T<:Real} = c * α
Base.:(//)(α::Cyclotomic, c::Real) = Cyclotomic(coeffs(α) .// c)
Base.:(/)(α::Cyclotomic, c::Real) = Cyclotomic(coeffs(α) ./ c)

function Base.div(α::Cyclotomic{T}, c::Number) where {T}
    RT = Base._return_type(div, (T, typeof(c)))
    rα = isone(conductor(α)) ? α : reduced_embedding(α)
    return div!(similar(rα, RT), rα, c)
end

###########################
# Complex arithmetic

function Base.promote_rule(
    ::Type{<:Cyclotomic{T}},
    ::Type{<:Complex{S}},
) where {T,S}
    return (TT = promote_type(T, S); Cyclotomic{TT,SparseVector{TT,Int}})
end

###########################
# Ring structure:

_enable_intermediate_reduction() = false

function _maybe_reduce(α::Cyclotomic{<:Rational{T}}) where {T<:Base.BitInteger}
    _enable_intermediate_reduction() || return α
    conductor(α) > 1 || return α

    k = (typemax(T) >> 4sizeof(T))
    for v in values(α)
        z = abs(v)
        if max(numerator(z), denominator(z)) > k
            return normalform!(reduced_embedding(α))
        end
    end
    return α
end

_maybe_reduce(α::Cyclotomic) = α

function common_embedding(α::Cyclotomic, β::Cyclotomic; reduced=true)
    if reduced
        α = reduced_embedding(α)
        β = reduced_embedding(β)
    end
    conductor(α) == conductor(β) && return α, β
    l = lcm(conductor(α), conductor(β))
    α = (l == conductor(α) ? α : embed(α, l))
    β = (l == conductor(β) ? β : embed(β, l))
    return α, β
end

function add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic)
    coeffs(out) .= coeffs(α) .+ coeffs(β)
    return out
end
function sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic)
    coeffs(out) .= coeffs(α) .- coeffs(β)
    return out
end

function mul!(out::Cyclotomic{T}, α::Cyclotomic, β::Cyclotomic) where {T}
    if out === α || out === β
        out = similar(out)
    end
    zero!(out)

    for (αe, αc) in exps_coeffs(α)
        for (βe, βc) in exps_coeffs(β)
            out[αe+βe] += αc * βc
        end
    end

    return out
end

for (op, fn) in ((:+, :add!), (:-, :sub!), (:*, :mul!))
    @eval begin
        function Base.$op(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
            α = _maybe_reduce(α)
            β = _maybe_reduce(β)
            α, β = common_embedding(α, β, reduced=false)
            @assert conductor(α) == conductor(β)
            U = promote_type(T, S)
            res = similar(α, U)
            try
                $fn(res, α, β)
            catch err
                @debug "overflow thrown, trying to reduce/normalize arguments!"
                if err isa OverflowError
                    α, β = normalform!.(common_embedding(α, β, reduced=true))
                    $fn(res, α, β)
                else
                    rethrow(err)
                end
            end
        end
    end
end

function Base.conj!(out::Cyclotomic, α::Cyclotomic, n::Integer = -1)
    zero!(out)
    for (e, c) in exps_coeffs(α)
        out[n*e] = c
    end
    return out
end

"""
    conj(α::Cyclotomic[, n::Integer=1])
Return the `n`-th conjugate of `α`, i.e. the image of `α` under the `n`-th
Frobenious homomorphism.

If `n` is co-prime to the conductor of `α` the map defines Galois automorphism.
Note that the default choice for `n=-1` corresponds to the standard complex
conjugation.
"""
function Base.conj(α::Cyclotomic, n::Integer = -1)
    return conj!(similar(α), α, n)
end

function galois_conj(α::Cyclotomic, n::Integer = -1)
    @assert gcd(n, conductor(α)) == 1
    return conj(α, n)
end

function Base.inv(α::Cyclotomic{T}) where {T}
    rα = isone(conductor(α)) ? α : reduced_embedding(α)
    RT = Base._return_type(inv, (T,))
    return inv!(similar(rα, RT), rα)
end

function inv!(
    out::Cyclotomic{T,<:AbstractSparseVector},
    α::Cyclotomic,
) where {T}
    copyto!(coeffs(out), coeffs(inv!(dense(zero!(out)), α)))
    return out
end

function inv!(
    out::Cyclotomic{T},
    α::Cyclotomic,
    tmp = zero(out),
    tmp2 = zero(out),
) where {T}
    if out === α
        out = one(out)
    else
        out = one!(out)
    end

    α = _maybe_reduce(α)

    let α = α
        basis_fb = zumbroich_viacomplement(conductor(α))

        # begin
        #     A1 = [conj(α, i) for i in _coprimes(conductor(α))]
        #     A2 = [normalform!(c, tmp, basis_forbidden = basis_fb) for c in A1]
        #     conjugates = unique!(coeffs, A2)
        #     normalform!(α, tmp, basis_forbidden=basis_fb)
        # end
        #
        # for c in conjugates
        #     c == α && continue
        #     mul!(tmp, out, c)
        #     copyto!(coeffs(out), coeffs(tmp))
        #     normalform!(out, tmp, basis_forbidden=basis_fb)
        # end

        for i in _coprimes(conductor(α))
            i < 2 && continue
            mul!(tmp2, out, conj!(tmp, α, i))
            copyto!(coeffs(out), coeffs(tmp2))
            normalform!(out, tmp, basis_forbidden = basis_fb)
        end

        # The idea:
        # out is the product of non-trivial Galois conjugates of α:
        # out = Π_{σ(Gal(𝕂(ζ_n)/𝕂)), σ≠id} σ(α)
        # since Π_{σ(Gal(𝕂(ζ_n)/𝕂))} σ(α) = norm_𝕂(α) ∈ 𝕂 we have
        # 1 = α·out/(α·out) = α · out/norm_𝕂(α), hence
        # α¯¹ = out/norm_𝕂(α)

        # however we don't necessarily take product of all of them as visible in the loop above;

        # @info reduced_embedding(out) == out

        norm_𝕂 = reduced_embedding(normalform!(out)*α)[0]
        w = out * inv(norm_𝕂)
        return _maybe_reduce(w)
    end

    return _maybe_reduce(out)
end
_coprimes(n::Integer) = BitSet(i for i in 2:n if gcd(i, n) == 1)

Base.:/(α::Cyclotomic, β::Cyclotomic) = α * inv(β)
