# Cyclotomics.jl

Cyclotomics package implements cyclotomic numbers which are sums of roots of unity.
The coefficients of the sum are in general taken from a ring.
E.g. the imaginary unit is represented by `E(4)`, the fourth root of `1`,
while algebraic number `(1 + √5)/2` can be written exactly as `E(5) + E(5)^4`.

In summary the package implements

* Cyclotomic numbers as structs based on `SparseVector`s,
* basic arithmetic on those: module and ring structures that take advantage of (lazy) normalization,
* a few predicates (e.g. `isreal`) and conversions to `float`/`Rational`/`Complex` numbers,
* Zumbroich basis (by three different methods), thread-safe and memoized.

## Example uses
```julia
julia> using Cyclotomics

julia> e = E(45) # 45-th root of unity
 ζ₄₅

julia> isone(E(5)^5) # 5-th root of unity to power 5 gives the unit
 true
```

### Normal forms
Consider the following element
```julia
julia> w = e + e^2 + e^8 + e^11 + e^17 + e^26 + e^29 + e^38 + e^44
 ζ₄₅ + ζ₄₅² + ζ₄₅⁸ + ζ₄₅¹¹ + ζ₄₅¹⁷ + ζ₄₅²⁶ + ζ₄₅²⁹ + ζ₄₅³⁸ + ζ₄₅⁴⁴
```

Since the vector space spanned by `45`-th roots of unity is of dimension less
than `44` not all roots are needed to express a cyclotomic number of degree `45`.
For example the following is a different way to write `w`:
```julia
julia> x = E(45) + E(45)^5 # or E(45) + E(9)
 ζ₄₅ + ζ₄₅² + ζ₄₅⁸ + ζ₄₅¹¹ + ζ₄₅¹⁷ + ζ₄₅²⁶ + ζ₄₅²⁹ + ζ₄₅³⁸ + ζ₄₅⁴⁴

julia> x == w
 true
```

And that's 9-th root of unity in its normal form (i.e. written in the canonical basis):
```julia
julia> E(45, 5) # == E(45)^5 == E(9)
 -ζ₉⁴-ζ₉⁷

```
### Computing with cyclotomics
We define module structures with different coefficients

```julia
julia> E(45, 5)
-ζ₉⁴-ζ₉⁷

julia> 3E(45, 5)
-3*ζ₉⁴ -3*ζ₉⁷

julia> 2.0E(45, 5) - E(9)
-1.0*ζ₉⁴ -1.0*ζ₉⁷

```
as well as conversions to standard julia types

```julia
julia> complex(2.0x)
3.5126250237210965 + 1.5639214212932075im

julia> float(E(3))
ERROR: InexactError: float(Float64,  1*E(3)^1)
Stacktrace:
[1] float(#unused#::Type{Float64}, α::Cyclotomic{Int64, SparseArrays.SparseVector{Int64, Int64}})
  @ Cyclotomics ~/.julia/dev/Cyclotomics/src/cycl.jl:168
[2] float(α::Cyclotomic{Int64, SparseArrays.SparseVector{Int64, Int64}})
  @ Cyclotomics ~/.julia/dev/Cyclotomics/src/cycl.jl:171
[3] top-level scope
  @ REPL[15]:1

julia> complex(E(3))
-0.4999999999999998 + 0.8660254037844387im

julia> float(E(3) + E(3)^2)
-1.0

julia> Rational(E(3) + E(3)^2)
-1//1

```

When possible we try to promote to Cyclotomics
```julia
julia> E(5) + im
-ζ₂₀ + ζ₂₀⁴-ζ₂₀⁹-ζ₂₀¹³-ζ₂₀¹⁷

julia> (1.0+2im) + E(5)
-2.0*ζ₂₀ -1.0*ζ₂₀⁸ -2.0*ζ₂₀⁹ -1.0*ζ₂₀¹² -2.0*ζ₂₀¹³ -1.0*ζ₂₀¹⁶ -2.0*ζ₂₀¹⁷

julia> (1.0+2.0im) - 2E(4)
1.0

julia> typeof(ans)
Cyclotomic{Float64, SparseArrays.SparseVector{Float64, Int64}}

julia> isreal((1.0+2.0im) - 2E(4))
true

```

However cyclotomic numbers can store non-rational algebraic numbers:

```julia
julia> z = E(5)^2+E(5)^3
 ζ₅² + ζ₅³

julia> isreal(z)
 true

julia> Rational(z)
ERROR: InexactError: Rational( 1*E(5)^2 + 1*E(5)^3)
Stacktrace:
 [1] Rational{Int64}(α::Cyclotomic{Int64, SparseArrays.SparseVector{Int64, Int64}})
   @ Cyclotomics ~/.julia/dev/Cyclotomics/src/cycl.jl:192
 [2] Rational(α::Cyclotomic{Int64, SparseArrays.SparseVector{Int64, Int64}})
   @ Cyclotomics ~/.julia/dev/Cyclotomics/src/cycl.jl:195
 [3] top-level scope
   @ none:1

julia> z ≈ (-1-sqrt(5))/2
true

```

### Low level constructors
One can also construct `Cyclotomic` directly from a vector, which is then used
as the underlying vector of coefficients. Here these are dense, while by default
`Cyclotomic` uses sparse storage.

```julia
julia> Cyclotomic(5, [0,1,0,0,0])
 ζ₅

julia> Cyclotomic(5, [0,0,1,0,0])
 ζ₅²

julia> Cyclotomic(5, [1,1,1,1,1]) # the sum of all roots == 0
 0

```


## `Cyclotomic`s: constructors and low level access

```@docs
Cyclotomic
Cyclotomics.E
Cyclotomics.conductor

Cyclotomics.exponents
Cyclotomics.coeffs
Base.values
Base.conj

Cyclotomics.dense
Cyclotomics.sparse
```

## Embeddings and normal forms

```@docs
Cyclotomics.reduced_embedding
Cyclotomics.embed

Cyclotomics.isnormalized
Cyclotomics.normalform!
Cyclotomics.normalform
```

### Internals

```@docs
Base.hash
```

## Technicalities: Zumbroich basis

One can naively represent cyclotomic number as a vector of `n` coefficients, corresponding to `n`-th root of identity. However since there are relations among them (e.g. the sum of all is equal to `0`), the actual dimension of the vector space is usually much smaller than `n`. Zumbroich basis is the set of `n`-th roots of unity which are linearly independent as vectors in the subspace (i.e. allow to express any cyclotomic number as sum of them).

There are three implementations provided by the package:
* [`zumbroich_plain`](https://github.com/kalmarek/Cyclotomics.jl/blob/76ceeb8822b1d63af2dab328c165a385b2af463a/src/zumbroich.jl#L14) following the description in the documentation of [GAP system](https://www.gap-system.org/Manuals/doc/ref/chap60_mj.html#X7F52BEA0862E06F2)
* [`zumbroich_direct`](https://github.com/kalmarek/Cyclotomics.jl/blob/76ceeb8822b1d63af2dab328c165a385b2af463a/src/zumbroich.jl#L38) which attempts to compute the basis directly
* [`zumbroich_viacomplement`](https://github.com/kalmarek/Cyclotomics.jl/blob/76ceeb8822b1d63af2dab328c165a385b2af463a/src/zumbroich.jl#L115) which computes the complement of the basis (the default).

My understanding and the implementation are based on the wonderful documenting comments at the top of [cyclotom.c](https://github.com/gap-system/gap/blob/master/src/cyclotom.c) from GAP project. Check them out!

!!! note The package uses function [`zumbroich_basis`](https://github.com/kalmarek/Cyclotomics.jl/blob/76ceeb8822b1d63af2dab328c165a385b2af463a/src/zumbroich.jl#L152) (which defaults to `zumbroich_viacomplement`) in its source code. To avoid recomputation the basis over and over the function is memoized for `Int` arguments.
