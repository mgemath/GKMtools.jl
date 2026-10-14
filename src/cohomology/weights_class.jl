@doc raw"""
    weight_class(G, e::Edge)

Return the equivariant first Chern class represented by the weight of `e`.
The result is an element of the localized equivariant coefficient ring of `G`.

# Examples
```jldoctest weight_class
julia> H7 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 7))
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (-7, -1, 0, 1)
4 -> 1 => (0, -1, -7, 1)
4 -> 3 => (1, 0, -1, 0)
Algorithmic connection for GKM graph with 4 nodes and valency 2

julia> weight_class(H7, Edge(3, 2))
-7*t1 - t2 + t4
```
"""
function weight_class(G, e::Edge)
  return _weight_class(G, e, gens_coeffRing(G))
end

function _weight_class(G, e::Edge, t::Vector{T}) where {T}
  res = zero(t[1])
  w = weight(G, e)
  for i in 1:rank_torus(G)
    res += w[i] * t[i]
  end
  return res
end

function _flag_weight_class(G, v::Int, i::Int, t::Vector{T}) where {T}
  res = zero(t[1])
  w = weight(G, v, i)
  for j in 1:rank_torus(G)
    res += w[j] * t[j]
  end
  return res
end

@doc raw"""
    unit_cohomology_ring(G)

Return the multiplicative identity in the localized GKM cohomology of `G`.
Its fixed-point restriction equals `1` at every vertex.
"""
function unit_cohomology_ring(G)
  return one(first(gens_cohomRing(G)))
end
