module Cyclotomics

import Primes
using SparseArrays

export coeffs, conductor, E, Cyclotomic

include("zumbroich.jl")
include("cycl.jl")
include("predicates.jl")
include("arithmetic.jl")
include("normalform.jl")
include("io.jl")

end # of module Cyclotomics
