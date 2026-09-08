"""
    quantum_schubert_matrix(Q, u; at_q1=false)
    quantum_schubert_matrix(products::AbstractVector; at_q1=false, p=nothing)

Return the matrix of quantum multiplication by `sigma_u`. Column `v`
contains `sigma_u ⋆ sigma_v`: entry `(w,v)` is the sum of all coefficients
`c*q^d` with output vertex `w`.

The first overload computes the complete product vector. The second accepts
a vector loaded from [`serialize_quantum_schubert_products`](@ref), without
needing a context or graph. The matrix size is `length(products)`.
Supply `p` when loading modular residues: the result is over GF(p), or its
Novikov polynomial ring. The context overload uses `Q.prime` automatically.
The saved product vector contains no prime metadata.

By default the result is an Oscar matrix over `QQ[q1,...,qr]`, with variables
in the order of the degree tuples. Set `at_q1=true` to sum all Novikov degrees
and return a matrix over `QQ`. An empty or entirely zero vector returns a
zero matrix over the selected coefficient field, since no Novikov coordinates
can be inferred.

# Examples
```julia
M = quantum_schubert_matrix(Q, 10)
M1 = quantum_schubert_matrix(Q, 10; at_q1=true)

using Serialization
products = deserialize("products-u10.jls")
M = quantum_schubert_matrix(products)
q = gens(base_ring(M))
M1 = quantum_schubert_matrix(products; at_q1=true)
```
"""
function quantum_schubert_matrix(products::AbstractVector; at_q1::Bool=false, p=nothing)
  K = _qs_field(p)
  n = length(products)
  rank_q = nothing
  for product in products, ((w,d),c) in product
    w isa Integer && 1 <= w <= n || throw(ArgumentError("output vertex must lie in 1:length(products)"))
    d isa Tuple && all(x -> x isa Integer && x >= 0, d) ||
      throw(ArgumentError("Novikov degree must be a tuple of nonnegative integers"))
    if isnothing(rank_q)
      rank_q = length(d)
    elseif length(d) != rank_q
      throw(ArgumentError("inconsistent numbers of Novikov coordinates"))
    end
  end
  r = something(rank_q, 0)
  if at_q1 || r == 0
    M = zero_matrix(K,n,n)
    for (v,product) in enumerate(products), ((w,d),c) in product
      M[w,v] += K(c)
    end
    return M
  end
  S,q = polynomial_ring(K, ["q$i" for i in 1:r])
  M = zero_matrix(S,n,n)
  for (v,product) in enumerate(products), ((w,d),c) in product
    term = S(c)
    for i in 1:r
      term *= q[i]^d[i]
    end
    M[w,v] += term
  end
  return M
end

function quantum_schubert_matrix(Q::QuantumSchubertContext,u; at_q1::Bool=false)
  return quantum_schubert_matrix(_qs_products(Q,u); at_q1, p=Q.prime)
end
