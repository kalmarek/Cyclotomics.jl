import Memoize
import LRUCache

Memoize.@memoize LRUCache.LRU{Tuple{Int},Tuple{BitSet,ForbiddenResidues{Int}}}(
    maxsize = 10000,
) zumbroich_viacomplement(n::Int) = zumbroich_viacomplement(n, Primes.factor(n))

Memoize.@memoize LRUCache.LRU{Tuple{Int},BitSet}(maxsize = 10000) function _coprimes(
    n::Int,
)
    return BitSet(i for i in 2:n if gcd(i, n) == 1)
end
