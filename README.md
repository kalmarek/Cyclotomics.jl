# Cyclotomics

![CI](https://github.com/kalmarek/Cyclotomics.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/kalmarek/Cyclotomics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/Cyclotomics.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kalmarek.github.io/Cyclotomics.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kalmarek.github.io/Cyclotomics.jl/dev)

## Example use:

```
using Cyclotomics

julia> e = E(45) # 45-th root of unity
 ζ₄₅

julia> 3e
  3*ζ₄₅

julia> w = e + e^2 + e^8 + e^11 + e^17 + e^26 + e^29 + e^38 + e^44
 ζ₄₅ + ζ₄₅² + ζ₄₅⁸ + ζ₄₅¹¹ + ζ₄₅¹⁷ + ζ₄₅²⁶ + ζ₄₅²⁹ + ζ₄₅³⁸ + ζ₄₅⁴⁴

julia> x = E(45) + E(45, 5) # this is equal to w
 ζ₄₅ + ζ₄₅² + ζ₄₅⁸ + ζ₄₅¹¹ + ζ₄₅¹⁷ + ζ₄₅²⁶ + ζ₄₅²⁹ + ζ₄₅³⁸ + ζ₄₅⁴⁴

julia> x == w
 true

julia> E(45, 5) # and that's 9-th root of unity in normal form
 -ζ₉⁴-ζ₉⁷

julia> 2.0*x
2.0*ζ₄₅ + 2.0*ζ₄₅² + 2.0*ζ₄₅⁸ + 2.0*ζ₄₅¹¹ + 2.0*ζ₄₅¹⁷ + 2.0*ζ₄₅²⁶ + 2.0*ζ₄₅²⁹ + 2.0*ζ₄₅³⁸ + 2.0*ζ₄₅⁴⁴

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

julia> isone(E(5)^5)
true

julia> E(5)^2+E(5)^3
 ζ₅² + ζ₅³

julia> (E(5)^2+E(5)^3) ≈ (-1-sqrt(5))/2
true

julia> Cyclotomic(5, [0,1,0,0,0])
 ζ₅

julia> Cyclotomic(5, [0,0,1,0,0])
 ζ₅²

julia> Cyclotomic(5, [1,1,1,1,1]) # the sum of all roots == 0
 0

```
