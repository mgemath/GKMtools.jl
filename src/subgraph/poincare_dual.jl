"""
    poincare_dual(S::AbstractGKMSubgraph)

Return the equivariant Poincaré dual cohomology class of the subgraph `S`
inside its ambient GKM graph.

For each vertex of the ambient graph that belongs to `S`, the class is formed
by multiplying the corresponding ambient cohomology generator by the product
of the weight classes of all ambient flags at that vertex that are not used by
the subgraph.

# Examples
```jldoctest poincare_dual
julia> P2 = projective_space(GKMGraph, 2);

julia> P1inP2 = subgraph_from_vertices(P2, [1, 2])
GKM subgraph of:
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, -1, 0)
3 -> 1 => (1, 0, -1)
3 -> 2 => (0, 1, -1)
Birkhoff-Grothendieck connection for GKM graph with 3 nodes and valency 2
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (1, -1, 0)
Algorithmic connection for GKM graph with 2 nodes and valency 1

julia> poincare_dual(P1inP2)
GKM class with restrictions:
[-t1 + t3, -t2 + t3, 0]

julia> is_gkm_spline(ans)
true
```
"""
function poincare_dual(S::AbstractGKMSubgraph)
  ambient = ambient_graph(S)
  H = get_cohomology(ambient)
  t = gens_coeffRing(ambient)

  result = zero(first(gens_cohomRing(ambient)))

  for local_vertex in 1:num_vertices(subgraph(S))
    ambient_vertex = vertex_to_ambient(S, local_vertex)
    contribution = unit_cohomology_ring(ambient)

    for ambient_flag in 1:length(flags(ambient, ambient_vertex))
      has_ambient_flag(S, local_vertex, ambient_flag) && continue
      contribution *= _flag_weight_class(ambient, ambient_vertex, ambient_flag, t)
    end

    result += contribution * gens_cohomRing(ambient)[ambient_vertex]
  end

  return result
end
