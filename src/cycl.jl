###############################################################################
#
#   Cyclotomics
abstract type AbstractCyclotomic <: Number end
"""
    Cyclotomic(coeffs::AbstractVector)
Element of `length(coeffs)`-th cyclotomic field with coefficients stored as `coeffs`.

To access the internals of a cyclotomic use API functions:
 * `conductor` - the conductor of a cyclotomic, i.e. the `n` used currently for storage. This might not be the minimal embeding field of a cyclotomic.
 * `coeffs` - returns a reference to the underlying vector containing the coefficients.
 * `getindex`/`setindex!` - use `α[i]` to access the coefficient at `i`-th power of the root of unity (in a circular fashion).
 * `values`/`exponents` - paired iterators over _non zero_ coefficients/exponents corresponding to _non-zero_ coefficients
 * `normalform!` - bring a cyclotomic to its unique representation as given by Zumbroich basis (also available in non-modifying form).

Iteration over non-zero coefficients in `Cyclotomic` is provided by `exps_coeffs`
which produces pairs `(exp, coeff)` of exponent and corresponding coefficient.
"""
struct Cyclotomic{T,A<:AbstractVector{T}} <: AbstractCyclotomic
    coeffs::A
end

Cyclotomic{T}(α::Cyclotomic) where {T} = Cyclotomic(convert.(T, coeffs(α)))

function Cyclotomic{T,A}(c::Cyclotomic{S,B}) where {T,S,A<:AbstractVector{T},B}
    return Cyclotomic{T,A}(coeffs(c))
end

Cyclotomic(c::Complex{T}) where {T} = Cyclotomic{T,SparseVector{T,Int}}(c)

Cyclotomic{T,A}(c::Real) where {T,A<:AbstractVector{T}} = Cyclotomic{T,A}(T[c])

function Cyclotomic{T,A}(c::Complex) where {T,A<:AbstractVector{T}}
    return Cyclotomic{T,A}(real(c)) + E(4) * Cyclotomic{T,A}(imag(c))
end

"""
    E(n[, i=1])
Return the `i`-th power of `n`-th root of unity with sparse vector as storage.
"""
function E(n, i = 1)
    k = Primes.totient(n)
    i = (0 <= i < n ? i : mod(i, n))
    coeffs = sparsevec([i + 1], [1], n)
    sizehint!(coeffs.nzind, k)
    sizehint!(coeffs.nzval, k)
    return Cyclotomic(coeffs)
end

####
#   Low-level interface

"""
    conductor(α::Cyclotomic)
Return the conductor, i.e. the degree of cyclotomic field `α` belongs to.
"""
conductor(α::Cyclotomic) = length(coeffs(α))

"""
    coeffs(α::Cyclotomic)
Return the coefficients of `α` as stored.
"""
coeffs(α::Cyclotomic) = α.coeffs

function _to_index(α::Cyclotomic, idx::Integer)
    # return mod(idx, conductor(α)) + 1
    0 <= idx < conductor(α) && return idx + 1
    conductor(α) <= idx && return (idx % conductor(α)) + 1
    return (idx % conductor(α)) + conductor(α) + 1
end

function Base.getindex(α::Cyclotomic, exp::Integer)
    return @inbounds α.coeffs[_to_index(α, exp)]
end

Base.getindex(α::Cyclotomic, itr) = [α[i] for i in itr]

function Base.setindex!(α::Cyclotomic, val, exp::Integer)
    @inbounds α.coeffs[_to_index(α, exp)] = val
    return val
end

struct ExpCoeffItr{C<:Cyclotomic}
    α::C
end

# general definitions for iteration
function Base.iterate(itr::ExpCoeffItr, state = 0)
    idx = findnext(!iszero, coeffs(itr.α), state + 1)
    idx === nothing && return nothing
    return (idx - 1, coeffs(itr.α)[idx]), idx
end

Base.IteratorSize(::Type{<:ExpCoeffItr}) = Base.HasLength()
Base.length(itr::ExpCoeffItr) = count(!iszero, coeffs(itr.α))
Base.eltype(::Type{<:ExpCoeffItr{<:Cyclotomic{T}}}) where {T} = Tuple{Int,T}

exps_coeffs(α::Cyclotomic) = ExpCoeffItr(α)
"""
    exponents(α::Cyclotomic)
Return an iterator over non-zero exponents of `α`, beginning at `0`-th one.
Matched iterator over coefficients is provided by @ref(values).
"""
exponents(α::Cyclotomic) = (e for (e, c) in exps_coeffs(α))

"""
    values(α::Cyclotomic)
Return an iterator over non-zero coefficients of `α`, beginning at `0`-th one.
Matched iterator over exponents is provided by @ref(exponents).
"""
Base.values(α::Cyclotomic) = (c for (e, c) in exps_coeffs(α))
Base.valtype(::Cyclotomic{T}) where {T} = T

function Base.similar(α::Cyclotomic, T::Type = valtype(α))
    return Cyclotomic(similar(coeffs(α), T))
end

function Base.similar(α::Cyclotomic, T::Type, n::Integer)
    return Cyclotomic(similar(coeffs(α), T, n))
end

"""
    dense(α::Cyclotomic)
Return a copy of `α` with coefficients stored in dense `Vector`.
"""
function dense(α::Cyclotomic{T}) where {T}
    return Cyclotomic{T,Vector{T}}(coeffs(α))
end

"""
    sparse(α::Cyclotomic)
Return a copy of `α` with coefficients stored in `SparseVector`.
"""
SparseArrays.sparse(α::Cyclotomic) = Cyclotomic(sparse(coeffs(α)))

function SparseArrays.droptol!(
    α::Cyclotomic{T,<:AbstractSparseVector},
    ε = eps(T) * 1e3,
) where {T}
    normalform!(α)
    coeffs(α) .= droptol!(coeffs(α), ε)
    return α
end

function roundcoeffs!(
    α::Cyclotomic{<:AbstractFloat},
    r::Base.RoundingMode = Base.RoundNearest;
    kwargs...,
)
    for i in eachindex(coeffs(α))
        coeffs(α)[i] = round(coeffs(α)[i], r; kwargs...)
    end
    return α
end

function Base.Complex{T}(α::Cyclotomic) where {T<:AbstractFloat}
    isone(conductor(α)) && return Complex{T}(α[0])
    rα = reduced_embedding(α)
    isone(conductor(rα)) && return Complex{T}(rα[0])
    n = conductor(rα)
    z = zero(Complex{T})
    for (e, c) in exps_coeffs(rα)
        γ = 2 * T(π) * T(e) / n
        z += T(c) * (cos(γ) + im * sin(γ))
    end
    return ifelse(rα == conj(rα), real(z) + zero(T)im, z)
end

Base.complex(α::Cyclotomic{T}) where {T} = Complex{float(T)}(α)
Base.complex(α::Cyclotomic{T}) where {T<:AbstractFloat} = Complex{T}(α)

function Base.float(::Type{T}, α::Cyclotomic) where {T<:AbstractFloat}
    αre, αim = real.(Complex{T}.(reim(α)))
    if abs(αim / αre) <= sqrt(eps(T)) || αim < eps(T)
        return αre
    end
    return throw(InexactError(:float, T, α))
end

Base.float(α::Cyclotomic) = float(float(valtype(α)), α)
Base.Float64(α::Cyclotomic) = float(Float64, α)

function _isreal(α::Cyclotomic)
    rα = reduced_embedding(α)
    if !isreal(rα) || any(!iszero, (rα[i] for i in 1:conductor(rα)-1))
        return (false, rα)
    else
        return (true, rα)
    end
end

function Base.Int(α::Cyclotomic)
    isone(conductor(α)) && return Int(α[0])
    flag, rα = _isreal(α)
    flag && return Int(rα[0])
    return throw(InexactError(:Int, Int, α))
end

function Base.Rational{T}(α::Cyclotomic) where {T}
    isone(conductor(α)) && return Rational{T}(α[0])
    flag, rα = _isreal(α)
    flag && return Rational{T}(rα[0])
    if isreal(rα)
        @error "The cyclotomic is real but it can not be converted to Rational: $rα ≈ $(float(rα))"
    end
    return throw(InexactError(:convert, Rational{T}, α))
end

Base.Rational(α::Cyclotomic{T}) where {T<:Integer} = Rational{T}(α)
Base.Rational(α::Cyclotomic{Rational{T}}) where {T} = Rational{T}(α)

Base.abs2(α::Cyclotomic) = α * conj(α)
Base.abs(α::Cyclotomic) = abs(complex(α))

Base.real(α::Cyclotomic) = (α + conj(α)) / 2
Base.imag(α::Cyclotomic) = -im * (α - conj(α)) / 2
