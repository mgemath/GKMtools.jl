### This file is dedicated to functions dealing with GKM graphs of different weight types
export convert_weights

@doc raw"""
    convert_weights(G::AbstractGKM_graph) -> AbstractGKM_graph{QQFieldElem}

It returns ``G``, but the character group will be embedded into a free ``\mathbb{Q}``-module.

# Examples
```jldoctest
julia> G_over_Z = generalized_gkm_flag(root_system(:A, 1))
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

julia> typeof(G_over_Z)
GKMtools.AbstractGKM_graph{ZZRingElem}

julia> G_over_Q = convert_weights(G_over_Z)
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

julia> typeof(G_over_Q)
GKMtools.AbstractGKM_graph{QQFieldElem}
```
"""
function convert_weights(G::AbstractGKM_graph)::AbstractGKM_graph{QQFieldElem}
  
  if typeof(G) == AbstractGKM_graph{QQFieldElem}
    return G
  end
  M = free_module(QQ, rank(G.M))
  rM = rank(M)
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()
  nv = n_vertices(G.g)
  weights_at_vertex = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}}(undef, nv)

  for k in keys(G.w)
    W[k] = sum(i -> QQ(G.w[k][i]) * gens(M)[i], 1:rM)
  end

  for i in 1:nv
    weights_at_vertex[i] = [sum(j -> QQ(w[j]) * gens(M)[j], 1:rM) for w in G.weights_at_vertex[i]]
  end

  return AbstractGKM_graph(G.g, G.labels, M, weights_at_vertex, G.edge_to_flag_index, G.flag_to_edge, W, G.equivariantCohomology, G.curveClasses, G.connection, G.QH_structure_consts, G.know_all_QH_structure_consts)
end

# Warning: this returns zero of the weight type.
# However, ZZ(0) == QQ(0) is true, so one must compare parent(ZZ(0)) vs parent(QQ(0)) instead.
function _get_weight_type(G::AbstractGKM_graph)::GKM_weight_type
    return zero(G.M)[1]
end