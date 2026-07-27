function generalized_gkm_schubert(::AbstractGKMGraph{R, V, F}, vertex_label; descending::Bool = true) where {V, R, F}
  error("It seems that the GKM graph is not of a generalized GKM flag variety.")
end



function generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_number::Int; descending::Bool = true) where {V <: AbstractFlagVertex, R, F}

  accept_vertex(v, w) = descending ? v < w : v > w

  accepted_vertices = Int[]

  vertices_G = labels(G)

  for i in 1:num_vertices(G)
    if accept_vertex(vertices_G[i], vertices_G[vertex_number])
      push!(accepted_vertices, i)
    end
  end

  push!(accepted_vertices, vertex_number)
  sort!(accepted_vertices)

  return subgraph_from_vertices(G, accepted_vertices)
end

function generalized_gkm_schubert(
  R::RootSystem,
  S::AbstractVector{<:RootSpaceElem},
  vertex_label::String,
  descending::Bool = true,
)
  G = generalized_gkm_flag(R, collect(S))
  return generalized_gkm_schubert(G, vertex_label; descending)
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, indices_of_S::Vector{RootSpaceElem}, vertex_label::String) -> AbstractGKM_subgraph

Let ``G`` be the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtools.generalized_gkm_flag).
This functions returns the subgraph of the variety ``G`` given by all Schubert cells corresponding to the points less or equal to the point labeled `vertex_label` in the Bruhat order.

# Examples
```jldoctest generalized_gkm_schubert
julia> R = root_system(:B, 2);

julia> generalized_gkm_schubert(R, "s1*s2")
GKM subgraph of:
GKM graph with 8 nodes, valency 4 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> s1 => (0, -1)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s2*s1 => (-1, 1)
s2 -> id => (0, -1)
s2 -> s2*s1 => (1, 1)
s1*s2 -> s1 => (-1, 0)
s1*s2 -> s1*s2*s1 => (1, 1)
s1*s2 -> s2 => (-1, 1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2*s1 => (-1, 0)
s2*s1*s2 -> s1*s2 => (0, -1)
s1*s2*s1*s2 -> s1 => (-1, -1)
s1*s2*s1*s2 -> s1*s2*s1 => (0, -1)
s1*s2*s1*s2 -> s2 => (-1, 0)
s1*s2*s1*s2 -> s2*s1*s2 => (-1, 1)
Birkhoff-Grothendieck connection for GKM graph with 8 nodes and valency 4
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (-1, 1)
s2 -> id => (0, -1)
s1*s2 -> s1 => (-1, 0)
s1*s2 -> s2 => (-1, 1)
Restricted connection for GKM graph with 4 nodes and valency 2
```

As before, the subset S can be a subset of simple roots or a subset of indices.
```jldoctest generalized_gkm_schubert
julia> generalized_gkm_schubert(R, [1], "s2")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s2 -> id => (0, -1)
s1*s2 -> id => (-1, 0)
s1*s2 -> s2 => (-1, 1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2 => (-1, 0)
s2*s1*s2 -> s1*s2 => (0, -1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s2 -> id => (0, -1)
Restricted connection for GKM graph with 2 nodes and valency 1

julia> S = simple_roots(R);

julia> generalized_gkm_schubert(R, [S[2]], "s1")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> id => (-1, -1)
s2*s1 -> s1 => (0, -1)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s1 => (-1, -1)
s1*s2*s1 -> s2*s1 => (-1, 1)
Cartan connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)
Restricted connection for GKM graph with 2 nodes and valency 1
```
"""
function generalized_gkm_schubert(
  R::RootSystem,
  indices_of_S::AbstractVector{<:Integer},
  vertex_label::String,
  descending::Bool = true,
)
  G = generalized_gkm_flag(R, indices_of_S)
  return generalized_gkm_schubert(G, vertex_label; descending)
end

function generalized_gkm_schubert(
  R::RootSystem,
  vertex_label::String,
  descending::Bool = true,
)
  return generalized_gkm_schubert(R, Int[], vertex_label; descending)
end

@doc raw"""
    generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_label::String; descending::Bool = true) where {V <: AbstractFlagVertex, R, F} -> GKMGraph

As above, but the generalized flag variety is already provided.
"""
function generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_label::String; descending::Bool = true) where {V <: AbstractFlagVertex, R, F}
  index = 0
  for v in 1:num_vertices(G)
    if label(G, v) == vertex_label
      index = v
      break
    end
  end
  @assert (index > 0) "Vertex $vertex_label not found"

  return generalized_gkm_schubert(G, index; descending)
end