@doc raw"""
    subgraph_from_vertices(G::AbstractGKMGraph, vertices_to_keep::AbstractVector{<:Integer})

Return the induced GKM subgraph of `G` on the specified vertices. Integer
indices are interpreted in the order supplied. Strings are matched against
`label(G, v)` and preserve the order of the requested labels.

# Examples
```jldoctest
julia> G = projective_space(GKM_graph, 3)
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> S = subgraph_from_vertices(G, [2, 3])
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
3 -> 2 => (0, -1, 1, 0)

julia> S.self
GKM graph with 2 nodes, valency 1 and axial function:
3 -> 2 => (0, -1, 1, 0)

julia> S.super
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> S2 = subgraph_from_vertices(G, [2, 3]; include_all_flags=true)
GKM graph with 2 nodes, valency 3 and axial function:
3 -> 2 => (0, -1, 1, 0)
Standalone flags:
2.1 => (-1, 1, 0, 0)
2.3 => (0, 1, 0, -1)
3.1 => (-1, 0, 1, 0)
3.3 => (0, 0, 1, -1)
```
"""
function subgraph_from_vertices(G::AbstractGKMGraph, vertices_to_keep::AbstractVector{<:Integer})
  vertex_map = Int.(vertices_to_keep)
  isempty(vertex_map) && throw(ArgumentError("the vertex array is empty"))
  all(v -> 1 <= v <= num_vertices(G), vertex_map) ||
    throw(ArgumentError("a vertex index lies outside the graph"))
  allunique(vertex_map) || throw(ArgumentError("vertex indices must be distinct"))
  return _induced_subgraph_from_vertices(G, vertex_map)
end

@doc raw"""
    subgraph_from_vertices(G::AbstractGKMGraph, labels_to_keep::AbstractVector{<:AbstractString})

As before, but the vertices are given by their labels.
"""
function subgraph_from_vertices(G::AbstractGKMGraph, labels_to_keep::AbstractVector{<:AbstractString})
  requested = String.(labels_to_keep)
  isempty(requested) && throw(ArgumentError("the label array is empty"))
  allunique(requested) || throw(ArgumentError("vertex labels must be distinct"))

  indices_by_label = Dict{String,Vector{Int}}()
  for v in vertices(G)
    push!(get!(indices_by_label, label(G, v), Int[]), v)
  end

  vertex_map = Vector{Int}(undef, length(requested))
  for (i, requested_label) in enumerate(requested)
    matches = get(indices_by_label, requested_label, Int[])
    isempty(matches) && throw(ArgumentError("unknown vertex label: $requested_label"))
    length(matches) == 1 || throw(ArgumentError(
      "vertex label is ambiguous: $requested_label",
    ))
    vertex_map[i] = only(matches)
  end
  return _induced_subgraph_from_vertices(G, vertex_map)
end

function _induced_subgraph_from_vertices(
  ambient::AbstractGKMGraph{R,V,F},
  vertex_map::Vector{Int},
) where {R,V,F}
  local_vertex = zeros(Int, num_vertices(ambient))
  for (local_index, global_vertex) in enumerate(vertex_map)
    local_vertex[global_vertex] = local_index
  end

  local_graph = Graph{Undirected}(length(vertex_map))
  local_flags = [F[] for _ in vertex_map]
  flag_map = [Int[] for _ in vertex_map]
  edge_flags = Dict{Edge,Tuple{Int,Int}}()

  for ambient_edge in edges(ambient)
    source = local_vertex[src(ambient_edge)]
    target = local_vertex[dst(ambient_edge)]
    (iszero(source) || iszero(target)) && continue

    ambient_source_flag, ambient_target_flag =
      _edge_flag_indices(ambient, ambient_edge)
    push!(local_flags[source], flags(ambient, src(ambient_edge))[ambient_source_flag])
    push!(local_flags[target], flags(ambient, dst(ambient_edge))[ambient_target_flag])
    push!(flag_map[source], ambient_source_flag)
    push!(flag_map[target], ambient_target_flag)
    add_edge!(local_graph, source, target)
    edge_flags[Edge(source, target)] = (
      length(local_flags[source]),
      length(local_flags[target]),
    )
  end

  data = GKMCombinatorialData{R,V,F}(
    local_graph,
    lattice(ambient),
    vertices_structure(ambient)[vertex_map],
    local_flags,
    edge_flags,
  )
  graph_of_subgraph = _graph_from_subgraph_data(
    ambient, data, vertex_map, flag_map,
  )
  return ambient isa OrbifoldGKMGraph ?
         OrbifoldGKMSubgraph(ambient, graph_of_subgraph, vertex_map, flag_map) :
         GKMSubgraph(ambient, graph_of_subgraph, vertex_map, flag_map)
end

function _subgraph_connection(data::GKMCombinatorialData{R}) where {R}
  isempty(data.edge_flags) && return empty_connection(R)
  return try
    # TODO: try to impelement a true restricted connection, when possible
    # build_gkm_connection(data; connection_type="Restricted")
    build_gkm_connection(data)
  catch error
    error isa ArgumentError || rethrow()
    empty_connection(R)
  end
end

function _graph_from_subgraph_data(
  ambient::GKMGraph{R,V,F},
  data::GKMCombinatorialData{R,V,F},
  vertex_map,
  flag_map,
) where {R,V,F}
  return GKMGraph{R,V,F}(
    data,
    _subgraph_connection(data),
    create_cohomology(rank(lattice(data)), length(vertex_map)),
    nothing,
  )
end

function _graph_from_subgraph_data(
  ambient::OrbifoldGKMGraph{R,V,F},
  data::GKMCombinatorialData{R,V,F},
  vertex_map,
  flag_map,
) where {R,V,F}
  vertex_isotropy = [
    let isotropy = ambient.vertex_isotropy[vertex_map[v]]
      OrbifoldVertexIsotropy(
        copy(isotropy.isotropy_group),
        isotropy.tangent_rep[:, flag_map[v]],
      )
    end
    for v in eachindex(vertex_map)
  ]
  flag_isotropy = [
    ambient.flag_isotropy[vertex_map[v]][flag_map[v]]
    for v in eachindex(vertex_map)
  ]
  return OrbifoldGKMGraph(
    data,
    vertex_isotropy,
    flag_isotropy,
    _subgraph_connection(data),
  )
end
