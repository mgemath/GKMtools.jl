@doc raw"""
    core(G::AbstractGKMGraph) -> GKMCombinatorialData

Return the combinatorial data underlying `G`.
"""
core(G::AbstractGKMGraph) = G.core

@doc raw"""
    graph(data::GKMCombinatorialData) -> Graph
    graph(G::AbstractGKMGraph) -> Graph

Return the undirected graph underlying combinatorial data or a GKM graph.
"""
graph(data::GKMCombinatorialData) = data.g
graph(G::AbstractGKMGraph) = graph(core(G))

@doc raw"""
    num_vertices(G::AbstractGKMGraph) -> Int

Return the number of vertices of `G`.
"""
num_vertices(G::AbstractGKMGraph) = nv(graph(G))

@doc raw"""
    num_edges(G::AbstractGKMGraph) -> Int

Return the number of edges of `G`.
"""
num_edges(G::AbstractGKMGraph) = ne(graph(G))

@doc raw"""
    lattice(data::GKMCombinatorialData)
    lattice(G::AbstractGKMGraph)

Return the character lattice of the torus acting on the GKM graph.
"""
lattice(data::GKMCombinatorialData) = data.M
lattice(G::AbstractGKMGraph) = lattice(core(G))

@doc raw"""
    labels(data::GKMCombinatorialData)
    labels(G::AbstractGKMGraph)

Return the vector of vertex-label objects, in vertex-index order. Use `label`
to obtain the display string of one vertex.
"""
labels(data::GKMCombinatorialData) = data.labels
labels(G::AbstractGKMGraph) = labels(core(G))

@doc raw"""
    flags(data::GKMCombinatorialData, v::Int)
    flags(G::AbstractGKMGraph, v::Int)

Return the ordered vector of flags at vertex `v`, including standalone flags.
"""
flags(data::GKMCombinatorialData, v::Int) = data.flags[v]
flags(G::AbstractGKMGraph, v::Int) = flags(core(G), v)

@doc raw"""
    edges(G::AbstractGKMGraph)

Return an iterator over the unoriented edges of `G`.
"""
edges(G::AbstractGKMGraph) = edges(graph(G))

@doc raw"""
    vertices(data::GKMCombinatorialData)
    vertices(G::AbstractGKMGraph)

Return an iterator over the vertex indices.
"""
vertices(core::GKMCombinatorialData) = vertices(graph(core))
vertices(G::AbstractGKMGraph) = vertices(core(G))

@doc raw"""
    label(G::AbstractGKMGraph, v::Int) -> String

Return the display label of vertex `v`.
"""
label(G::AbstractGKMGraph, v::Int) = get_string(labels(G)[v])

@doc raw"""
    degree(G::AbstractGKMGraph, v::Int) -> Int

Return the graph-theoretic degree of vertex `v`. Standalone flags do not
contribute to this degree.
"""
degree(G::AbstractGKMGraph, v::Int) = degree(graph(G), v)


@doc raw"""
    find_vertex_index(vertex_label::String, G::AbstractGKMGraph) -> Int

Return the index of the vertex whose display label is `vertex_label`. An
`AssertionError` is raised if no vertex has that label.
"""
function find_vertex_index(Vertexlabel::String, G::AbstractGKMGraph)
  index = 0
  for i in 1:num_vertices(G)
    if Vertexlabel == label(G, i)
      index = i
      break
    end
  end
  @assert (index > 0) "label not found"
  return index
end


@doc raw"""
    compact_flags(G::AbstractGKMGraph, v::Int) -> Vector{Int}

Return the sorted indices of the flags at vertex `v` which belong to edges.
Standalone flags are omitted.
"""
function compact_flags(G::AbstractGKMGraph, v::Int)
  
  if degree(G, v) == length(flags(G, v))
    return collect(1:length(flags(G, v)))
  end

  ans = Vector{Int}(undef, degree(G, v))
  index = 1
  for e in edges(G)
    if src(e) == v
      i, _ = core(G).edge_flags[e]
      ans[index] = i
      index += 1
    elseif dst(e) == v
      _, j = core(G).edge_flags[e]
      ans[index] = j
      index += 1
    end
  end
  
  return sort(ans)
end

@doc raw"""
    weight(data::GKMCombinatorialData, v::Int, i::Int)
    weight(G::AbstractGKMGraph, v::Int, i::Int)
    weight(data::GKMCombinatorialData, e::Edge)
    weight(G::AbstractGKMGraph, e::Edge)

Return an axial weight. The `(v, i)` forms return the weight of the `i`-th
flag at vertex `v`. The edge forms orient the weight from `src(e)` to `dst(e)`;
reversing `e` therefore negates the result.
"""
function weight(core::GKMCombinatorialData, v::Int, i::Int)
  return core.flags[v][i].weight
end

function weight(G::AbstractGKMGraph, v::Int, i::Int)
  return weight(core(G), v, i)
end

function weight(core::GKMCombinatorialData, e::Edge)
  _e = e
  sign = 1
  if !haskey(core.edge_flags, e)
    # Try the opposite edge for undirected graphs
    _e = Edge(dst(e), src(e))
    sign = -1
    
    if !haskey(core.edge_flags, _e)
      error("Edge $e not found in edge_flags")
    end
  end
  i, _ = core.edge_flags[_e]
  v = src(_e)
  return sign * weight(core, v, i)
end

function weight(G::AbstractGKMGraph, e::Edge)
  return weight(core(G), e)
end


@doc raw"""
    other_vertex(G::AbstractGKMGraph, v::Int, i::Int) -> Int

Return the vertex opposite to the `i`-th flag at vertex `v`.

This is the vertex on the other side of the corresponding edge. Return `0` if
the flag is standalone.
"""
function other_vertex(G::AbstractGKMGraph, v::Int, i::Int)
  for edge in edges(G)
    if src(edge) == v
      source_flag, _ = _edge_flag_indices(G, edge)
      source_flag == i && return dst(edge)
    elseif dst(edge) == v
      _, target_flag = _edge_flag_indices(G, edge)
      target_flag == i && return src(edge)
    end
  end
  # throw(ArgumentError("flag $i at vertex $v does not belong to an edge"))
  return 0
end

function _edge_flag_indices(G::AbstractGKMGraph, e::Edge)
  return _edge_flag_indices(core(G), e)
end

function _edge_flag_indices(data::GKMCombinatorialData, e::Edge)
  if haskey(data.edge_flags, e)
    return data.edge_flags[e]
  end
  target, source = data.edge_flags[reverse(e)]
  return source, target
end

@doc raw"""
    print_labels(G::AbstractGKMGraph)

Print all strings of the label of the vertices of `G` with their vertex numbers.

# Examples
```jldoctest
julia> R = root_system(:A, 2);

julia> G = generalized_gkm_flag(R)
GKM graph with 6 nodes, valency 3 and axial function:
s1 -> id => (-1, 1, 0)
s2*s1 -> s1 => (0, -1, 1)
s1*s2*s1 -> id => (-1, 0, 1)
s1*s2*s1 -> s2*s1 => (-1, 1, 0)
s2 -> id => (0, -1, 1)
s2 -> s2*s1 => (1, 0, -1)
s1*s2 -> s1 => (-1, 0, 1)
s1*s2 -> s1*s2*s1 => (0, 1, -1)
s1*s2 -> s2 => (-1, 1, 0)
Birkhoff-Grothendieck connection for GKM graph with 6 nodes and valency 3

julia> print_labels(G)
1 -> id
2 -> s1
3 -> s2*s1
4 -> s1*s2*s1
5 -> s2
6 -> s1*s2
```
"""
function print_labels(G::AbstractGKMGraph)
  for i in 1:num_vertices(G)
    println("$i -> ", label(G, i))
  end
end

@doc raw"""
    number_vertex_of_label(G::AbstractGKMGraph, label::String) -> Int

Return the vertex number of a given label. An `AssertionError` is raised if no
vertex has that label.
# Examples
```jldoctest
julia> R = root_system(:A, 2);

julia> G = generalized_gkm_flag(R);

julia> number_vertex_of_label(G, "id")
1

julia> number_vertex_of_label(G, "s1")
2
```
"""
function number_vertex_of_label(G::AbstractGKMGraph, label::String)
  for i in 1:num_vertices(G)
    if label == label(G, i)
      return i
    end
  end
  error("label $label not found")
end