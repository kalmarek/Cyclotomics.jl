# Cyclotomics.jl

The `Cyclotomics` package implements cyclotomic numbers which are sums of roots of unity. E.g. the imaginary unit is represented by `E(4)`, the fourths root of `1`. The coefficients of the sum are in general taken from a ring.

In summary the package implements

* Cyclotomic numbers as structs based on `AbstractVectors` with basic functionality,
* basic arithmetic: module and ring structures that takes advantage of (lazy) normalization,
* a few predicates and conversions to floats/rationals/complex numbers.
* Zumbroich basis (by three different methods), thread-safe and memoized,

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
* `zumbroich_plain` following the description in the documentation of [GAP system](https://www.gap-system.org/Manuals/doc/ref/chap60_mj.html#X7F52BEA0862E06F2)
* `zumbroich_direct` which attempts to compute the basis directly
* `zumbroich_viacomplement` which computes the complement of the basis (the default).

My understanding and the implementation are based on the wonderful documenting comments at the top of [cyclotom.c](https://github.com/gap-system/gap/blob/master/src/cyclotom.c) from GAP project. Check them out!

Note: The package uses function `zumbroich_basis` (which defaults to `zumbroich_viacomplement`) in its source code. To avoid recomputation the basis over and over the function is memoized for `Int`s.
