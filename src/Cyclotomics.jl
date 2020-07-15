module Cyclotomics

using Primes
using SparseArrays

export zumbroich, Cyclotomic, exponents, conductor, E

include("zumbroich.jl")
include("cycl.jl")
include("predicates.jl")
include("arithmetic.jl")
include("normalform.jl")
include("io.jl")

end # of module Cyclotomics