"""
    serialize_quantum_schubert_products(path, Q, u; show_bar=false)

Compute all ordinary quantum products `sigma_u ⋆ sigma_v` for fixed `u`,
save them with Julia's `Serialization.serialize`, and return the saved vector.

The vector entry `products[v]` is a sparse dictionary
`(output_vertex, degree_tuple) => BigInt`. Only this vector is stored:
no graph, context, cache, representatives, or other metadata. Coefficients
are native Julia integers, so loading does not require Oscar or GKMtools.
Indices and degree coordinates retain the conventions of `Q`.
For modular contexts, integers are residues in `0:Q.prime-1`; the prime is
not stored. Supply it when reconstructing a matrix from the saved vector.

Set `show_bar=true` to display progress over the second factor.
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
  path::AbstractString, Q::QuantumSchubertContext, u; show_bar::Bool=false,
)
  products = _qs_products(Q, u; show_bar)
  Serialization.serialize(path, products)
  return products
end

function _qs_products(Q::QuantumSchubertContext, u; show_bar::Bool=false)
  u = _qs_vertex(Q, u)
  products = Vector{Dict{Tuple{Int,_QSDegree},BigInt}}(undef,length(Q.representatives))
  progress = show_bar ? ProgressMeter.Progress(length(products)) : nothing
  for v in eachindex(products)
    products[v] = Dict{Tuple{Int,_QSDegree},BigInt}(
      key => _qs_integer(Q,c) for (key,c) in quantum_schubert_product(Q,u,v)
    )
    isnothing(progress) || ProgressMeter.next!(progress)
  end
  return products
end
