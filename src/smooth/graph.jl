###############################################################################
# CORE GKM GRAPH STRUCTURE (RECOMMENDED)
###############################################################################

struct GKMGraph{R, V, F} <: AbstractGKMGraph{R, V, F}
  # where R is the coefficient ring for the character lattice, and V is the type of vertex labels

  core::GKMCombinatorialData{R, V, F}

  ###########################################################################
  # CACHED GEOMETRIC DATA
  ###########################################################################

  # Connection (canonical or algorithmic)
  connection::Connection{R}

  # Equivariant cohomology
  cohomology::GKMCohomology

  # Curve classes (H₂)
  H2::GKM_H2

  # Quantum / GW data
  quantum::Union{Nothing,GKM_Quantum}

  function GKMGraph{R,V,F}(core::GKMCombinatorialData{R,V,F}, connection::Connection{R}, cohomology::GKMCohomology, quantum::Union{Nothing,GKM_Quantum}) where {R,V,F}
    H2 = _GKM_second_homology(core)
    return new{R,V,F}(core, connection, cohomology, H2, quantum)
  end

  function GKMGraph{R,V,F}(core::GKMCombinatorialData{R,V,F}, connection::Connection{R}, cohomology::GKMCohomology, H2::GKM_H2, quantum::Union{Nothing,GKM_Quantum}) where {R,V,F}
    return new{R,V,F}(core, connection, cohomology, H2, quantum)
  end

  function GKMGraph{R,V,F}(core::GKMCombinatorialData{R,V,F}, connection::Connection{R}, cohomology::GKMCohomology) where {R,V,F}
    H2 = _GKM_second_homology(core)
    return new{R,V,F}(core, connection, cohomology, H2, nothing)
  end
end

function Base.show(io::IO, G::GKMGraph)
  print(io, "GKMGraph with $(nv(G.core.g)) stacky vertices")
end

function Base.show(io::IO, ::MIME"text/plain", G::GKMGraph)
  print(
    io, "GKM graph with $(n_vertices(graph(G))) nodes, valency $(valency(G)) and axial function:"
  )
  for e in edges(G)
    print(io, "\n$(label(G, src(e))) -> $(label(G, dst(e))) => $(weight(G, e))")
  end

  if !is_compact(G)
    print(io, "\nStandalone flags:")
  end

  print_connection(io, G; verbose = false)

  # print standalone flags if any:
  # is_compact(G) && return nothing
  # print(io, "\nStandalone flags:")
  # for v in 1:n_vertices(G.g)
  #   for (i, w) in enumerate(G.weights_at_vertex[v])
  #     !isnothing(G.flag_to_edge[v][i]) && continue
  #     print(io, "\n$(label(G, v)).$i => $w")
  #   end
  # end
end

function gkm_graph(core::GKMCombinatorialData{R, V, F}) where {R, V, F}
  connection = build_gkm_connection(core)
  M = lattice(core)
  g = graph(core)
  cohomology = create_cohomology(rank(M), nv(g))
  return GKMGraph{R,V,F}(core, connection, cohomology, nothing)
end

# function GKMGraph(
#   g::Graph,
#   M::AbstractAlgebra.Generic.FreeModule{R},
#   labels::Vector{V},
#   flags::Vector{Vector{F}},
#   edge_flags::Dict{Edge,Tuple{Int,Int}}
# ) where {R, V, F}
#   core = GKMCombinatorialData{R, V, F}(
#     g,
#     M,
#     labels,
#     flags,
#     edge_flags
#   )
#   return GKMGraph{R, V, F}(
#     core,
#     nothing,   # connection
#     nothing,   # H2
#     nothing,   # cohomology
#     nothing,    # quantum
#   )
# end

# Interface

# function Base.show(io::IO, gkm::GKMGraph{R}) where R
#     print(io, "GKMGraph{$R} with $(nv(gkm.g)) vertices and $(ne(gkm.g)) edges")

#     # Indicate which lazy fields are computed
#     status = String[]
#     gkm.connection !== nothing && push!(status, "Connection")
#     gkm.cohomology !== nothing && push!(status, "Cohomology")
#     gkm.quantum    !== nothing && push!(status, "Quantum")

#     if !isempty(status)
#         print(io, " (Computed: ", join(status, ", "), ")")
#     end
# end

# @doc raw"""
#     gkm_graph(g, labels, M, w; check=true, checkLabels=true) -> GKMGraph

# Create a GKM graph from the given data.

# # Arguments
# - `g::Graph`: An undirected graph.
# - `labels::Vector{String}`: Labels for vertices.
# - `M::AbstractAlgebra.Generic.FreeModule{R}`: Character lattice.
# - `w::Dict{Edge, FreeModuleElem{R}}`: Axial function (one orientation per edge).

# # Notes
# - Only constructs the **combinatorial GKM graph**.
# - Connections, H₂, cohomology etc. are NOT computed here.
# # Example
# Let us construct the GKM graph of the projective line. First of all, we create a graph with two vertices, and one edge.
# ```jldoctest first_GKM_graph
# julia> g = Graph{Undirected}(2)
# Undirected graph with 2 nodes and no edges

# julia> add_edge!(g, 1, 2);
# ```
# Let us define our array of labels.

# ```jldoctest first_GKM_graph
# julia> labels = ["a", "b"];
# ```
# Now, we create the character group. We take a free module of rank 2 over the integers.

# ```jldoctest first_GKM_graph
# julia> M = free_module(ZZ, 2)
# Free module of rank 2 over ZZ
# ```
# We create the axial function. It is a dictionary from the set of edges to the character group. This time we have only one edge.

# ```jldoctest first_GKM_graph
# julia> e = first(edges(g));

# julia> w = Dict(e => gens(M)[1] - gens(M)[2])
# Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}} with 1 entry:
#   Edge(2, 1) => (1, -1)
# ```
# Finally, we create the GKM graph.
# ```jldoctest first_GKM_graph
# julia> gkm_graph(g, labels, M, w)
# GKM graph with 2 nodes, valency 1 and axial function:
# b -> a => (1, -1)
# ```

# !!! warning
#     1. Do not change the number of verices after this.
#     2. Don't modify the underlying OSCAR graph directly after this. Use the functions of this package instead.
#     3. All edges should be added immediately after calling this function and not changed afterwards.

# !!! note
#     After you have added all edges using `add_edge!`, you may use `initialize!` to calculate the GKM connection (if it is unique) and the curve classes. If you don't do this, those data will be calculated whenever required for the first time.

# """
# function gkm_graph(
#   g::Graph,
#   labels::Vector{String},
#   M::AbstractAlgebra.Generic.FreeModule{R},
#   w::Dict{Edge,AbstractAlgebra.Generic.FreeModuleElem{R}};
#   check::Bool=true,
#   checkLabels::Bool=true,
# ) where {R}
#   n_vert = nv(g)

#   ###########################################################################
#   # 1. VALIDATION
#   ###########################################################################
#   if check
#     @req n_vert ≥ 1 "GKM graph needs at least one vertex"
#     @req length(labels) == n_vert "Mismatch in number of labels"
#     @req allunique(labels) "Labels must be unique"
#     @req Set(edges(g)) == keys(w) "Axial function not defined on all edges"
#     @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"

#     # regular valency
#     deg = degree(g, 1)
#     @req all(v -> degree(g, v) == deg, 1:n_vert) "Graph must have constant valency"
#   end

#   if checkLabels
#     @req all(l -> !occursin(r"[>\[\]]", l), labels) "Forbidden characters in labels"
#   end

#   ###########################################################################
#   # 2. COMPLETE AXIAL FUNCTION (ADD REVERSE EDGES)
#   ###########################################################################
#   w_full = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}()

#   for e in edges(g)
#     @req haskey(w, e) "Missing weight for edge $e"
#     w_full[e] = w[e]
#     w_full[reverse(e)] = -w[e]
#   end

#   ###########################################################################
#   # 3. BUILD FLAGS (CORE DESIGN CHANGE)
#   ###########################################################################
#   flags = Vector{Vector{GKMFlag{R}}}(undef, n_vert)
#   edge_flags = Dict{Edge,Tuple{Int,Int}}()

#   for v in 1:n_vert
#     nbrs = neighbors(g, v)
#     deg = length(nbrs)

#     flags[v] = Vector{GKMFlag{R}}(undef, deg)

#     for (i, u) in enumerate(nbrs)
#       e = Edge(v, u)
#       flags[v][i] = GKMFlag(v, w_full[e], e)
#     end
#   end

#   ###########################################################################
#   # 4. BUILD EDGE ↔ FLAG CORRESPONDENCE
#   ###########################################################################
#   for e in edges(g)
#     v, u = src(e), dst(e)

#     i = findfirst(f -> f.edge == e, flags[v])
#     j = findfirst(f -> f.edge == reverse(e), flags[u])

#     @req i !== nothing && j !== nothing "Inconsistent flag construction"

#     edge_flags[e] = (i, j)
#   end

#   ###########################################################################
#   # 5. CONSTRUCT GRAPH (NO DERIVED DATA!)
#   ###########################################################################
#   return GKMGraph(g, M, labels, flags, edge_flags, nothing, nothing, nothing, nothing)
# end