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

function Base.div(α::Cyclotomic, c::Number)
    T = typeof(div(α[0], c))
    return div!(similar(α, T), normalform!(α), c)
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

_enable_intermediate_normalization() = false

function _maybe_normalize!(
    α::Cyclotomic{<:Rational{T}},
) where {T<:Base.BitInteger}

    # if false
    if _enable_intermediate_normalization() && !isnormalized(α)
        k = (typemax(T) >> 4 * sizeof(T))
        for v in values(α)
            z = abs(v)
            if max(numerator(z), denominator(z)) > k
                normalform!(α)
                break
            end
        end
    end
    return α
end

_maybe_normalize!(α::Cyclotomic) = α

function add!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic)
    coeffs(out) .= coeffs(α) .+ coeffs(β)
    return _maybe_normalize!(out)
end
function sub!(out::Cyclotomic, α::Cyclotomic, β::Cyclotomic)
    coeffs(out) .= coeffs(α) .- coeffs(β)
    return _maybe_normalize!(out)
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

    return _maybe_normalize!(out)
end

for (op, fn) in ((:+, :add!), (:-, :sub!), (:*, :mul!))
    @eval begin
        function Base.$op(α::Cyclotomic{T}, β::Cyclotomic{S}) where {T,S}
            if conductor(α) == conductor(β)
                return $fn(similar(α, promote_type(T, S)), α, β)
            else
                l = lcm(conductor(α), conductor(β))
                return $op(embed(α, l), embed(β, l))
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
    rα = reduced_embedding(α)
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
    out::Cyclotomic{T,<:DenseVector},
    α::Cyclotomic,
    tmp = similar(out),
    tmp2 = similar(out),
) where {T}
    if out === α
        out = one(out)
    else
        out = one!(out)
    end

    ilead = inv(maximum(abs, coeffs(α)))
    T <: AbstractFloat &&
        ilead < eps(T) &&
        @warn "Invering element with large lead: $(maximum(abs, coeffs(α)))" α
    let α = α * ilead
        basis, fb = zumbroich_viacomplement(conductor(α))
        lb = length(basis)
        conjugates_counter = 0

        for i in 2:conductor(α)-1
            conjugates_counter == lb - 1 && break
            any(x -> gcd(i, first(x)) > 1, fb) && continue
            conjugates_counter += 1
            mul!(tmp2, out, conj!(tmp, α, i))
            copyto!(coeffs(out), coeffs(tmp2))
        end

        # out is now the product of non-trivial Galois conjugates of α:
        # out = Π_{σ(Gal(𝕂(ζ_n)/𝕂)), σ≠id} σ(α)
        # since Π_{σ(Gal(𝕂(ζ_n)/𝕂))} σ(α) = norm_𝕂(α) ∈ 𝕂 we have
        # 1 = α·out/(α·out) = α · out/norm_𝕂(α), hence
        # α¯¹ = out/norm_𝕂(α)

        norm_𝕂 = reduced_embedding(mul!(tmp2, out, α))
        # norm_𝕂 should be real by now

        if T <: AbstractFloat
            float(imag(norm_𝕂)) <= sqrt(eps(T)) * conductor(α) ||
                @warn "norm_𝕂  should be real, but it has imaginary part of magnitude $(float(imag(norm_𝕂)))"
            norm_α = float(real(norm_𝕂))
            # @info α norm_α
            # @show float.(reim(norm_𝕂))
            out = mul!(out, out, inv(norm_α))
        else
            @assert conductor(norm_𝕂) == 1 "$norm_ℚ" conductor(norm_𝕂)
            norm_α = norm_𝕂[0]
            out = mul!(out, out, inv(norm_α))
        end
    end

    mul!(out, out, ilead)
    return _maybe_normalize!(out)
end

Base.:/(α::Cyclotomic, β::Cyclotomic) = α * inv(β)
