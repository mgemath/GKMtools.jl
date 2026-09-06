"""
    serialize_quantum_schubert_products(path, Q, u)

Compute all ordinary quantum products `sigma_u ⋆ sigma_v` for fixed `u`,
save them with Julia's `Serialization.serialize`, and return the saved vector.

The vector entry `products[v]` is a sparse dictionary
`(output_vertex, degree_tuple) => BigInt`. Only this vector is stored:
no graph, context, cache, representatives, or other metadata. Coefficients
are native Julia integers, so loading does not require Oscar or GKMtools.
Indices and degree coordinates retain the conventions of `Q`.

All possible Novikov degrees are included. Computation finishes before the
destination is opened; an existing file is overwritten.

# Example
```julia
serialize_quantum_schubert_products("products-u3.jls", Q, 3)

using Serialization
products = deserialize("products-u3.jls")
products[2] # sigma_3 ⋆ sigma_2
```
"""
function serialize_quantum_schubert_products(
  path::AbstractString, Q::QuantumSchubertContext, u,
)
  products = _qs_products(Q, u)
  Serialization.serialize(path, products)
  return products
end

function _qs_products(Q::QuantumSchubertContext, u)
  u = _qs_vertex(Q, u)
  return [
    Dict{Tuple{Int,_QSDegree},BigInt}(
      key => BigInt(ZZ(c)) for (key,c) in quantum_schubert_product(Q,u,v)
    )
    for v in eachindex(Q.representatives)
  ]
end
