###############################################################################
# CORE GKM GRAPH STRUCTURE (RECOMMENDED)
###############################################################################

mutable struct GKMGraph{R} <: AbstractGKMGraph{R}

  ###########################################################################
  # CORE COMBINATORIAL DATA
  ###########################################################################

  g::Graph

  # Character lattice (typically ZZ or QQ module)
  M::AbstractAlgebra.Generic.FreeModule{R}

  # Vertex labels
  labels::Vector{String}

  # flags[v] = list of flags at vertex v
  flags::Vector{Vector{GKMFlag{R}}}

  # edge -> (flag index at source, flag index at target)
  edge_flags::Dict{Edge,Tuple{Int,Int}}

  ###########################################################################
  # LAZY CACHED GEOMETRIC DATA
  ###########################################################################

  # Connection (canonical or algorithmic)
  connection::Union{Nothing,GKMConnection}

  # Curve classes (H₂)
  H2::Union{Nothing,GKM_H2}

  # Equivariant cohomology
  cohomology::Union{Nothing,GKM_Cohomology}

  # Quantum / GW data
  quantum::Union{Nothing,GKM_Quantum}
end

function GKMGraph(
  g::Graph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  labels::Vector{String},
  flags::Vector{Vector{GKMFlag{R}}},
  edge_flags::Dict{Edge,Tuple{Int,Int}},
) where {R}
  return GKMGraph{R}(
    g,
    M,
    labels,
    flags,
    edge_flags,
    nothing,   # connection
    nothing,   # H2
    nothing,   # cohomology
    nothing,    # quantum
  )
end

# Interface

graph(g::GKMGraph) = g.g
lattice(g::GKMGraph) = g.M
flags(g::GKMGraph, v::Int) = g.flags[v]
edge_flags(g::GKMGraph) = g.edge_flags

num_vertices(g::GKMGraph) = nv(g.g)
num_edges(g::GKMGraph) = ne(g.g)

function weight(g::GKMGraph{R}, e::Edge) where {R}
  i, _ = g.edge_flags[e]
  v = src(e)
  return g.flags[v][i].weight
end

@doc raw"""
    gkm_graph(g, labels, M, w; check=true, checkLabels=true) -> GKMGraph

Create a GKM graph from the given data.

# Arguments
- `g::Graph`: An undirected graph.
- `labels::Vector{String}`: Labels for vertices.
- `M::AbstractAlgebra.Generic.FreeModule{R}`: Character lattice.
- `w::Dict{Edge, FreeModuleElem{R}}`: Axial function (one orientation per edge).

# Notes
- Only constructs the **combinatorial GKM graph**.
- Connections, H₂, cohomology etc. are NOT computed here.
# Example
Let us construct the GKM graph of the projective line. First of all, we create a graph with two vertices, and one edge.
```jldoctest first_GKM_graph
julia> g = Graph{Undirected}(2)
Undirected graph with 2 nodes and no edges

julia> add_edge!(g, 1, 2);
```
Let us define our array of labels.

```jldoctest first_GKM_graph
julia> labels = ["a", "b"];
```
Now, we create the character group. We take a free module of rank 2 over the integers.

```jldoctest first_GKM_graph
julia> M = free_module(ZZ, 2)
Free module of rank 2 over ZZ
```
We create the axial function. It is a dictionary from the set of edges to the character group. This time we have only one edge.

```jldoctest first_GKM_graph
julia> e = first(edges(g));

julia> w = Dict(e => gens(M)[1] - gens(M)[2])
Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}} with 1 entry:
  Edge(2, 1) => (1, -1)
```
Finally, we create the GKM graph.
```jldoctest first_GKM_graph
julia> gkm_graph(g, labels, M, w)
GKM graph with 2 nodes, valency 1 and axial function:
b -> a => (1, -1)
```

!!! warning
    1. Do not change the number of verices after this.
    2. Don't modify the underlying OSCAR graph directly after this. Use the functions of this package instead.
    3. All edges should be added immediately after calling this function and not changed afterwards.

!!! note
    After you have added all edges using `add_edge!`, you may use `initialize!` to calculate the GKM connection (if it is unique) and the curve classes. If you don't do this, those data will be calculated whenever required for the first time.

"""
function gkm_graph(
  g::Graph,
  labels::Vector{String},
  M::AbstractAlgebra.Generic.FreeModule{R},
  w::Dict{Edge,AbstractAlgebra.Generic.FreeModuleElem{R}};
  check::Bool=true,
  checkLabels::Bool=true,
) where {R}
  n_vert = nv(g)

  ###########################################################################
  # 1. VALIDATION
  ###########################################################################
  if check
    @req n_vert ≥ 1 "GKM graph needs at least one vertex"
    @req length(labels) == n_vert "Mismatch in number of labels"
    @req allunique(labels) "Labels must be unique"
    @req Set(edges(g)) == keys(w) "Axial function not defined on all edges"
    @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"

    # regular valency
    deg = degree(g, 1)
    @req all(v -> degree(g, v) == deg, 1:n_vert) "Graph must have constant valency"
  end

  if checkLabels
    @req all(l -> !occursin(r"[>\[\]]", l), labels) "Forbidden characters in labels"
  end

  ###########################################################################
  # 2. COMPLETE AXIAL FUNCTION (ADD REVERSE EDGES)
  ###########################################################################
  w_full = Dict{Edge,typeof(first(values(w)))}()

  for e in edges(g)
    @req haskey(w, e) "Missing weight for edge $e"
    w_full[e] = w[e]
    w_full[reverse(e)] = -w[e]
  end

  ###########################################################################
  # 3. BUILD FLAGS (CORE DESIGN CHANGE)
  ###########################################################################
  flags = Vector{Vector{GKMFlag{R}}}(undef, n_vert)
  edge_flags = Dict{Edge,Tuple{Int,Int}}()

  for v in 1:n_vert
    nbrs = neighbors(g, v)
    deg = length(nbrs)

    flags[v] = Vector{GKMFlag{R}}(undef, deg)

    for (i, u) in enumerate(nbrs)
      e = Edge(v, u)
      flags[v][i] = GKMFlag(v, w_full[e], e)
    end
  end

  ###########################################################################
  # 4. BUILD EDGE ↔ FLAG CORRESPONDENCE
  ###########################################################################
  for e in edges(g)
    v, u = src(e), dst(e)

    i = findfirst(f -> f.edge == e, flags[v])
    j = findfirst(f -> f.edge == reverse(e), flags[u])

    @req i !== nothing && j !== nothing "Inconsistent flag construction"

    edge_flags[e] = (i, j)
  end

  ###########################################################################
  # 5. CONSTRUCT GRAPH (NO DERIVED DATA!)
  ###########################################################################
  return GKMGraph(g, M, labels, flags, edge_flags, nothing, nothing, nothing, nothing)
end