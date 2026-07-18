function generalized_gkm_schubert(::AbstractGKMGraph{R, V, F}, vertex_label; descending::Bool = true) where {V, R, F}
  error("It seems that the GKM graph is not of a generalized GKM flag variety.")
end

function generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_label::String, descending::Bool = true) where {V <: AbstractFlagVertex, R, F}
  index = 0
  for v in 1:num_vertices(G)
    if label(G, v) == vertex_label
      index = v
      break
    end
  end
  @assert (index > 0) "Vertex $vertex_label not found"

  return generalized_gkm_schubert(G, index, descending)
end

function generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_number::Int, descending::Bool = true) where {V <: AbstractFlagVertex, R, F}

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
  return generalized_gkm_schubert(G, vertex_label, descending)
end

function generalized_gkm_schubert(
  R::RootSystem,
  indices_of_S::AbstractVector{<:Integer},
  vertex_label::String,
  descending::Bool = true,
)
  G = generalized_gkm_flag(R, indices_of_S)
  return generalized_gkm_schubert(G, vertex_label, descending)
end

function generalized_gkm_schubert(
  R::RootSystem,
  vertex_label::String,
  descending::Bool = true,
)
  return generalized_gkm_schubert(R, Int[], vertex_label, descending)
end
