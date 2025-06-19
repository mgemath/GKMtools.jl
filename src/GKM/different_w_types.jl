### This file is dedicated to functions dealing with GKM graphs of different weight types
@doc raw"""
    convert_weights(G::GKM_graph) -> GKM_graph{QQFieldElem}

It returns ``G``, but the character group will be embedded into a free ``\mathbb{Q}``-module.

# Examples
```jldoctest
julia> G_over_Z = generalized_gkm_flag(root_system(:A, 1))
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

julia> typeof(G_over_Z)
GKM_graph{ZZRingElem}

julia> G_over_Q = convert_weights(G_over_Z)
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

julia> typeof(G_over_Q)
GKM_graph{QQFieldElem}
```
"""
function convert_weights(G::GKM_graph)::GKM_graph{QQFieldElem}
  
  if typeof(G) == GKM_graph{QQFieldElem}
    return G
  end
  M = free_module(QQ, rank(G.M))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()

  for k in keys(G.w)
    W[k] = sum(i -> QQ(G.w[k][i]) * gens(M)[i], 1:rank(M))
  end

  return GKM_graph(G.g, G.labels, M, W, G.equivariantCohomology, G.H2, G.connection, G.QH_structure_consts, G.know_all_QH_structure_consts)  
end

function _get_weight_type(G::GKM_graph)::Weight
    return zero(G.M)[1]
end