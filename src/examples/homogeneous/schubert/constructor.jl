function generalized_gkm_schubert(::AbstractGKMGraph{R, V, F}, vertex_label; descending::Bool = true) where {V, R, F}
  error("It seems that the GKM graph is not of a generalized GKM flag variety.")
end



function generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_number::Int; descending::Bool = true) where {V <: AbstractFlagVertex, R, F}

  accept_vertex(v, w) = descending ? v < w : v > w

  accepted_vertices = Int[]

  vertices_G = vertices_structure(G)

  for i in 1:num_vertices(G)
    if accept_vertex(vertices_G[i], vertices_G[vertex_number])
      push!(accepted_vertices, i)
    end
  end

  push!(accepted_vertices, vertex_number)
  sort!(accepted_vertices)

  # Check regularity before the smooth GKM homology construction, whose
  # assumptions need not hold for a singular Schubert moment graph.
  accepted = Set(accepted_vertices)
  degrees = [count(in(accepted), neighbors(graph(G), v)) for v in accepted_vertices]
  schubert_label = label(G, vertex_number)
  all(==(first(degrees)), degrees) || throw(ArgumentError(
    "The Schubert variety indexed by $schubert_label is not rationally smooth: " *
    "its moment graph does not have constant valency. " *
    "generalized_gkm_schubert supports smooth Schubert varieties only.",
  ))

  return try
    subgraph_from_vertices(G, accepted_vertices)
  catch err
    err isa ArgumentError || rethrow()
    throw(ArgumentError(
      "Cannot construct the smooth GKM graph of the Schubert variety indexed by $schubert_label. " *
      "Its moment graph has constant valency, but rational smoothness does not imply smoothness. " *
      "Singular Schubert varieties are not supported by generalized_gkm_schubert. " *
      "Underlying error: $(sprint(showerror, err))",
    ))
  end
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
    generalized_gkm_schubert(R::RootSystem, S::AbstractVector, vertex_label::String, descending::Bool = true) -> AbstractGKMSubgraph
    generalized_gkm_schubert(R::RootSystem, vertex_label::String, descending::Bool = true) -> AbstractGKMSubgraph

Let ``G`` be the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtools.generalized_gkm_flag).
Return the induced subgraph of ``G`` on the points less than or equal to
`vertex_label` in the Bruhat order. Omitting `S` uses the full flag variety.
With `descending=false`, use the points greater than or equal to `vertex_label`
instead, corresponding to the opposite Schubert variety.

!!! warning "Smoothness requirement"
    This constructor supports smooth Schubert varieties. If the moment graph
    has nonconstant vertex valency, the variety is not rationally smooth and
    an `ArgumentError` is thrown before constructing the GKM subgraph.
    Constant valency (rational smoothness) does not guarantee smoothness.
    For a rationally smooth but singular variety, the smooth GKM construction
    can still fail, for example because Chern weights are not proportional to
    an edge weight. Such `ArgumentError`s include this distinction and the
    underlying failure. A successful return is not a smoothness test;
    smoothness must be established separately before using the geometric
    invariants of the returned graph.

# Examples
```jldoctest generalized_gkm_schubert
julia> R = root_system(:B, 2);

julia> generalized_gkm_schubert(R, "s1*s2")
GKM subgraph of:
GKM graph with 8 nodes, valency 4 and axial function:
s1 -> id => (1, -1)
s2*s1 -> s1 => (0, 1)
s1*s2*s1 -> id => (1, 0)
s1*s2*s1 -> s2*s1 => (1, -1)
s2 -> id => (0, 1)
s2 -> s2*s1 => (-1, -1)
s1*s2 -> s1 => (1, 0)
s1*s2 -> s1*s2*s1 => (-1, -1)
s1*s2 -> s2 => (1, -1)
s2*s1*s2 -> id => (1, 1)
s2*s1*s2 -> s2*s1 => (1, 0)
s2*s1*s2 -> s1*s2 => (0, 1)
s1*s2*s1*s2 -> s1 => (1, 1)
s1*s2*s1*s2 -> s1*s2*s1 => (0, 1)
s1*s2*s1*s2 -> s2 => (1, 0)
s1*s2*s1*s2 -> s2*s1*s2 => (1, -1)
Birkhoff-Grothendieck connection for GKM graph with 8 nodes and valency 4
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (1, -1)
s2 -> id => (0, 1)
s1*s2 -> s1 => (1, 0)
s1*s2 -> s2 => (1, -1)
Algorithmic connection for GKM graph with 4 nodes and valency 2
```

As before, the subset S can be a subset of simple roots or a subset of indices.
```jldoctest generalized_gkm_schubert
julia> generalized_gkm_schubert(R, [1], "s2")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s2 -> id => (0, 1)
s1*s2 -> id => (1, 0)
s1*s2 -> s2 => (1, -1)
s2*s1*s2 -> id => (1, 1)
s2*s1*s2 -> s2 => (1, 0)
s2*s1*s2 -> s1*s2 => (0, 1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s2 -> id => (0, 1)
Algorithmic connection for GKM graph with 2 nodes and valency 1

julia> S = simple_roots(R);

julia> generalized_gkm_schubert(R, [S[2]], "s1")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s1 -> id => (1, -1)
s2*s1 -> id => (1, 1)
s2*s1 -> s1 => (0, 1)
s1*s2*s1 -> id => (1, 0)
s1*s2*s1 -> s1 => (1, 1)
s1*s2*s1 -> s2*s1 => (1, -1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (1, -1)
Algorithmic connection for GKM graph with 2 nodes and valency 1
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
  return generalized_gkm_schubert(R, Int[], vertex_label, descending)
end

@doc raw"""
    generalized_gkm_schubert(G::AbstractGKMGraph{R, V, F}, vertex_label::String; descending::Bool = true) where {V <: AbstractFlagVertex, R, F} -> AbstractGKMSubgraph
    generalized_gkm_schubert(G::AbstractGKMGraph, vertex_number::Int; descending::Bool = true) -> AbstractGKMSubgraph

As above, but the generalized flag variety is already provided. The same
smoothness requirement and `ArgumentError` diagnostics apply to both vertex
labels and vertex numbers, including when `descending=false`.

# Examples
```jldoctest generalized_gkm_schubert
julia> R = root_system(:B, 2);

julia> generalized_gkm_schubert(R, "s1*s2")
GKM subgraph of:
GKM graph with 8 nodes, valency 4 and axial function:
s1 -> id => (1, -1)
s2*s1 -> s1 => (0, 1)
s1*s2*s1 -> id => (1, 0)
s1*s2*s1 -> s2*s1 => (1, -1)
s2 -> id => (0, 1)
s2 -> s2*s1 => (-1, -1)
s1*s2 -> s1 => (1, 0)
s1*s2 -> s1*s2*s1 => (-1, -1)
s1*s2 -> s2 => (1, -1)
s2*s1*s2 -> id => (1, 1)
s2*s1*s2 -> s2*s1 => (1, 0)
s2*s1*s2 -> s1*s2 => (0, 1)
s1*s2*s1*s2 -> s1 => (1, 1)
s1*s2*s1*s2 -> s1*s2*s1 => (0, 1)
s1*s2*s1*s2 -> s2 => (1, 0)
s1*s2*s1*s2 -> s2*s1*s2 => (1, -1)
Birkhoff-Grothendieck connection for GKM graph with 8 nodes and valency 4
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (1, -1)
s2 -> id => (0, 1)
s1*s2 -> s1 => (1, 0)
s1*s2 -> s2 => (1, -1)
Algorithmic connection for GKM graph with 4 nodes and valency 2
```

```jldoctest generalized_gkm_schubert
julia> generalized_gkm_schubert(R, [1], "s2")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s2 -> id => (0, 1)
s1*s2 -> id => (1, 0)
s1*s2 -> s2 => (1, -1)
s2*s1*s2 -> id => (1, 1)
s2*s1*s2 -> s2 => (1, 0)
s2*s1*s2 -> s1*s2 => (0, 1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s2 -> id => (0, 1)
Algorithmic connection for GKM graph with 2 nodes and valency 1

julia> S = simple_roots(R);

julia> generalized_gkm_schubert(R, [S[2]], "s1")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s1 -> id => (1, -1)
s2*s1 -> id => (1, 1)
s2*s1 -> s1 => (0, 1)
s1*s2*s1 -> id => (1, 0)
s1*s2*s1 -> s1 => (1, 1)
s1*s2*s1 -> s2*s1 => (1, -1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (1, -1)
Algorithmic connection for GKM graph with 2 nodes and valency 1

```
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
