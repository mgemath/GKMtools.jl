abstract type AbstractGKMSubgraph{R,V,F} end

"""
    GKMSubgraph(ambient, subgraph, vertex_map, flag_map)

A smooth or orbifold GKM subgraph embedded in an ambient graph.

`vertex_map[v]` is the ambient vertex corresponding to subgraph vertex `v`.
`flag_map[v][i]` is the corresponding ambient flag.

Forward and dense reverse maps are stored, giving constant-time translation
between subgraph and ambient vertices and flags. A zero in a reverse map means
that the ambient object is not part of the subgraph.
"""
struct GKMSubgraph{
  R,
  V,
  F,
  A<:AbstractGKMGraph{R,V,F},
  S<:AbstractGKMGraph{R,V,F},
} <: AbstractGKMSubgraph{R,V,F}
  ambient::A
  graph::S

  vertex_to_ambient::Vector{Int}
  ambient_to_vertex::Vector{Int}

  flag_to_ambient::Vector{Vector{Int}}
  ambient_to_flag::Vector{Vector{Int}}
end

function GKMSubgraph(
  ambient::A,
  sub::S,
  vertex_to_ambient::AbstractVector{<:Integer},
  flag_to_ambient::AbstractVector{
    <:AbstractVector{<:Integer}
  },
) where {
  R,
  V,
  F,
  A<:AbstractGKMGraph{R,V,F},
  S<:AbstractGKMGraph{R,V,F},
}
  ambient_is_orbifold =
    ambient isa AbstractOrbifoldGKMGraph

  subgraph_is_orbifold =
    sub isa AbstractOrbifoldGKMGraph

  ambient_is_orbifold == subgraph_is_orbifold ||
    throw(
      ArgumentError(
        "smooth and orbifold GKM graphs cannot be mixed",
      ),
    )

  n_local = num_vertices(sub)
  n_ambient = num_vertices(ambient)

  length(vertex_to_ambient) == n_local ||
    throw(
      DimensionMismatch(
        "expected $n_local vertex images, got " *
        "$(length(vertex_to_ambient))",
      ),
    )

  length(flag_to_ambient) == n_local ||
    throw(
      DimensionMismatch(
        "expected flag maps for $n_local vertices, got " *
        "$(length(flag_to_ambient))",
      ),
    )

  vertex_map = Int.(vertex_to_ambient)

  all(v -> 1 <= v <= n_ambient, vertex_map) ||
    throw(
      ArgumentError(
        "a vertex image lies outside the ambient graph",
      ),
    )

  allunique(vertex_map) ||
    throw(
      ArgumentError("the vertex map is not injective"),
    )

  reverse_vertices = zeros(Int, n_ambient)

  for (local_vertex, ambient_vertex) in enumerate(vertex_map)
    reverse_vertices[ambient_vertex] = local_vertex
  end

  flag_map = Vector{Vector{Int}}(undef, n_local)
  reverse_flags = Vector{Vector{Int}}(undef, n_local)

  for local_vertex in 1:n_local
    ambient_vertex = vertex_map[local_vertex]
    local_flags = Int.(flag_to_ambient[local_vertex])

    n_local_flags =
      length(flags(sub, local_vertex))

    n_ambient_flags =
      length(flags(ambient, ambient_vertex))

    length(local_flags) == n_local_flags ||
      throw(
        DimensionMismatch(
          "expected $n_local_flags flag images at vertex " *
          "$local_vertex, got $(length(local_flags))",
        ),
      )

    all(i -> 1 <= i <= n_ambient_flags, local_flags) ||
      throw(
        ArgumentError(
          "a flag image at vertex $local_vertex lies " *
          "outside the ambient flag set",
        ),
      )

    allunique(local_flags) ||
      throw(
        ArgumentError(
          "the flag map at vertex $local_vertex is not injective",
        ),
      )

    flag_map[local_vertex] = local_flags
    reverse_flags[local_vertex] =
      zeros(Int, n_ambient_flags)

    for (local_flag, ambient_flag) in enumerate(local_flags)
      reverse_flags[local_vertex][ambient_flag] =
        local_flag
    end
  end

  # Validate edge incidence once, so translating an edge later requires
  # only two array lookups.
  for e in edges(sub)
    ambient_edge = Edge(
      vertex_map[src(e)],
      vertex_map[dst(e)],
    )

    has_edge(graph(ambient), ambient_edge) ||
      throw(
        ArgumentError(
          "the image of subgraph edge $e is not an ambient edge",
        ),
      )

    local_source_flag, local_target_flag =
      _edge_flag_indices(sub, e)

    ambient_source_flag, ambient_target_flag =
      _edge_flag_indices(ambient, ambient_edge)

    flag_map[src(e)][local_source_flag] ==
      ambient_source_flag ||
      throw(
        ArgumentError(
          "the source flag of edge $e has an inconsistent image",
        ),
      )

    flag_map[dst(e)][local_target_flag] ==
      ambient_target_flag ||
      throw(
        ArgumentError(
          "the target flag of edge $e has an inconsistent image",
        ),
      )
  end

  return GKMSubgraph{R,V,F,A,S}(
    ambient,
    sub,
    vertex_map,
    reverse_vertices,
    flag_map,
    reverse_flags,
  )
end

@doc raw"""
    ambient_graph(G::GKMSubgraph)

Return the ambient graph of `G`.
"""
ambient_graph(G::GKMSubgraph) = G.ambient

@doc raw"""
    subgraph(G::GKMSubgraph)

Return the subgraph of `G` as a smooth or orbifold GKM graph, thus forgetting the ambient space.
"""
subgraph(G::GKMSubgraph) = G.graph

vertex_to_ambient(
  G::GKMSubgraph,
  vertex::Int,
) = G.vertex_to_ambient[vertex]

ambient_to_vertex(
  G::GKMSubgraph,
  vertex::Int,
) = G.ambient_to_vertex[vertex]

flag_to_ambient(
  G::GKMSubgraph,
  vertex::Int,
  flag::Int,
) = G.flag_to_ambient[vertex][flag]

ambient_to_flag(
  G::GKMSubgraph,
  local_vertex::Int,
  ambient_flag::Int,
) = G.ambient_to_flag[local_vertex][ambient_flag]

function edge_to_ambient(
  G::GKMSubgraph,
  e::Edge,
)
  return Edge(
    G.vertex_to_ambient[src(e)],
    G.vertex_to_ambient[dst(e)],
  )
end

function edge_from_ambient(
  G::GKMSubgraph,
  e::Edge,
)
  source = G.ambient_to_vertex[src(e)]
  target = G.ambient_to_vertex[dst(e)]

  iszero(source) && return nothing
  iszero(target) && return nothing

  local_edge = Edge(source, target)

  return has_edge(graph(G.graph), local_edge) ?
         local_edge :
         nothing
end

function has_ambient_vertex(
  G::GKMSubgraph,
  vertex::Int,
)
  return !iszero(G.ambient_to_vertex[vertex])
end

function has_ambient_flag(
  G::GKMSubgraph,
  local_vertex::Int,
  ambient_flag::Int,
)
  return !iszero(
    G.ambient_to_flag[local_vertex][ambient_flag],
  )
end

function Base.show(io::IO, G::GKMSubgraph)

  kind =
    G.graph isa AbstractOrbifoldGKMGraph ?
    "Orbifold GKM" :
    "Smooth GKM"

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "$kind GKM subgraph")
  else
    # nested printing allowed, preferably terse
    print(io, "$kind GKM subgraph with $(n_vertices(subgraph(G))) nodes and valency $(valency(subgraph(G)))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::GKMSubgraph)

  println(io, "GKM subgraph of:")
  show(io, MIME"text/plain"(), ambient_graph(G))
  println(io, "\nSubgraph:")
  show(io, MIME"text/plain"(), subgraph(G))
end
