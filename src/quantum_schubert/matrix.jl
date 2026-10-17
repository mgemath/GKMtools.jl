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
```jldoctest
julia> R = root_system(:C, 3);

julia> S = [1, 3];

julia> Q = quantum_schubert_context(R, S);

julia> Q.representatives[2]
s2

julia> M = quantum_schubert_matrix(Q, 2; at_q1=true)
[0   0   0   0   0   0   1   0   0   0   0   0]
[1   0   0   0   0   0   0   0   1   0   0   0]
[0   1   0   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   0   0   0   1   0]
[0   0   1   1   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   1]
[0   0   0   0   1   2   0   0   0   0   0   0]
[0   0   0   0   2   1   0   0   0   0   0   0]
[0   0   0   0   0   0   1   1   0   0   0   0]
[0   0   0   0   0   0   0   1   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   0   0]
[0   0   0   0   0   0   0   0   0   0   1   0]

julia> chi = characteristic_polynomial(M)
x^12 - 26*x^7 - 27*x^2

julia> eigenvalues_with_multiplicities(algebraic_closure(QQ), M)
11-element Vector{Tuple{QQBarFieldElem, Int64}}:
 ({a5: 1.93318}, 1)
 ({a1: 0}, 2)
 ({a1: -1.00000}, 1)
 ({a4: 0.809017 + 0.587785*im}, 1)
 ({a4: 0.809017 - 0.587785*im}, 1)
 ({a5: 0.597386 + 1.83857*im}, 1)
 ({a5: 0.597386 - 1.83857*im}, 1)
 ({a4: -0.309017 + 0.951057*im}, 1)
 ({a4: -0.309017 - 0.951057*im}, 1)
 ({a5: -1.56398 + 1.13630*im}, 1)
 ({a5: -1.56398 - 1.13630*im}, 1)
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
