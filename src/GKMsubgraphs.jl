import Oscar.has_edge
import Oscar.has_vertex

@doc"""
    gkm_subgraph_from_vertices(gkm::AbstractGKM_graph, vertices::Vector{Int64}; include_all_flags::Bool=false) -> AbstractGKM_subgraph

Return the GKM subgraph induced by the given vertices.

# Arguments
- `gkm`: The supergraph
- `vertices`: Vector of vertex indices to include in the subgraph
- `include_all_flags`: If `true`, includes all flags at the subgraph vertices. If `false` (default), only includes flags corresponding to edges contained in the subgraph.

!!! note
    1. This does not check if the result is a valid GKM graph (use `isvalid` for that).
    2. If possible, the subgraph will be endowed with the connection induced from the supergraph.
    3. If `include_all_flags=true` then the valency of the subgraph will be identical with the valency of the given graph.

# Example
```jldoctest subgr_from_vert
julia> G = projective_space(GKM_graph, 3)
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> S = gkm_subgraph_from_vertices(G, [2, 3])
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

julia> S2 = gkm_subgraph_from_vertices(G, [2, 3]; include_all_flags=true).self
GKM graph with 2 nodes, valency 3 and axial function:
3 -> 2 => (0, -1, 1, 0)
Standalone flags:
2.1 => (-1, 1, 0, 0)
2.3 => (0, 1, 0, -1)
3.1 => (-1, 0, 1, 0)
3.3 => (0, 0, 1, -1)
```
"""
function gkm_subgraph_from_vertices(gkm::AbstractGKM_graph, vertices::Vector{Int64}; include_all_flags::Bool=false) :: AbstractGKM_subgraph

  @req all(v -> v>0, vertices) "Vertex index must be positive"
  nv = n_vertices(gkm.g)
  @req all(v -> v <= nv, vertices) "Vertex index exceeds number of vertices"

  return _gkm_subgraph_from_vertices(gkm, unique(sort(vertices)), include_all_flags)
end

function _gkm_subgraph_from_vertices(gkm::AbstractGKM_graph, vDict::Vector{Int64}, include_all_flags::Bool) :: AbstractGKM_subgraph

  # Build sub_flags: which flags from the supergraph to include at each vertex
  subnv = length(vDict)
  sub_flags = Vector{Vector{Int64}}(undef, subnv)

  if include_all_flags
    # Include all flags at each vertex
    val_super = valency(gkm)
    for i in 1:subnv
      sub_flags[i] = collect(1:val_super)
    end
  else
    # Only include flags corresponding to edges that connect vertices in the subgraph
    for i in 1:subnv
      v_super = vDict[i]
      sub_flags[i] = Int64[]

      # Check each flag at this vertex in the supergraph
      for flag_idx in 1:length(gkm.weights_at_vertex[v_super])
        edge_super = gkm.flag_to_edge[v_super][flag_idx]

        # Include if it's an edge connecting to another vertex in the subgraph
        if !isnothing(edge_super)
          # v_neighbor = (src(edge_super) == v_super) ? dst(edge_super) : src(edge_super)
          @req src(edge_super) == v_super "edge-flag relations broken. Check if the gkm graph is valid!"
          
          v_neighbor = dst(edge_super)
          if v_neighbor in vDict
            push!(sub_flags[i], flag_idx)
          end
        end
      end
    end
  end

  return gkm_subgraph_from_flags(gkm, vDict, sub_flags)
end

# If the supergraph's connection is compatible with the subgraph, infer it to the subgraph and return true.
# Else return false. False is also returned if the subgraph's connection is already set.
function _infer_GKM_connection!(gkmSub::AbstractGKM_subgraph)::Bool

  # Try to get connection, but handle errors (e.g., for non-uniform valency graphs)
  con = nothing
  try
    con = get_connection(gkmSub.super)
  catch e
    # If connection computation fails (e.g., due to non-uniform valency), just skip
    return false
  end

  if !isnothing(con) && is_compatible_with_connection(gkmSub, con; printDiagnostics=false)

    oldCon = con.con
    oldA = con.a
    newCon = Dict{Edge, Vector{Int64}}()
    newA = Dict{Edge, Vector{ZZRingElem}}()
    vDict = gkmSub.vDict

    for e in edges(gkmSub.self.g)
      # Number of flags at src(e) and dst(e) (may be different for non-uniform valency)
      val_src = length(gkmSub.flagDict[src(e)])
      val_dst = length(gkmSub.flagDict[dst(e)])

      newCon[e] = Vector{Int64}(undef, val_src)
      newCon[reverse(e)] = Vector{Int64}(undef, val_dst)
      newA[e] = Vector{ZZRingElem}(undef, val_src)
      newA[reverse(e)] = Vector{ZZRingElem}(undef, val_dst)

      eSup = Edge(vDict[src(e)], vDict[dst(e)])

      for i in 1:val_src
        # Get the flag index in the supergraph using flagDict
        iSup = gkmSub.flagDict[src(e)][i]

        # Get the connected flag in the supergraph
        jSup = oldCon[eSup][iSup]

        # Map back to subgraph flag index
        j = findfirst(k -> gkmSub.flagDict[dst(e)][k] == jSup, 1:val_dst)
        if isnothing(j)
          error("Connection incompatible: connected flag not in subgraph")
        end

        newCon[e][i] = j
        newA[e][i] = oldA[eSup][iSup]

        # Also set the reverse connection
        newCon[reverse(e)][j] = i
        newA[reverse(e)][j] = oldA[reverse(eSup)][jSup]
      end
    end

    newConObj = GKM_connection(gkmSub.self, newCon, newA)
    set_connection!(gkmSub.self, newConObj)
    return true
  end
  return false
end

@doc"""
    gkm_subgraph_from_vertices(gkm::AbstractGKM_graph, vertexLabels::Vector{String}; include_all_flags::Bool=false) -> AbstractGKM_subgraph

As before, but the vertices are given by their labels.
"""
function gkm_subgraph_from_vertices(gkm::AbstractGKM_graph, vertexLabels::Vector{String}; include_all_flags::Bool=false) :: AbstractGKM_subgraph

  @req all(l -> l in gkm.labels, vertexLabels) "Label not found"

  vertices::Vector{Int64} = indexin(vertexLabels, gkm.labels) # need to specify Vector{Int64} as indexin returns vector of Union{Nothing, Int64}.
  return gkm_subgraph_from_vertices(gkm, vertices; include_all_flags=include_all_flags)
end

@doc"""
    gkm_subgraph_from_flags(gkm::AbstractGKM_graph, sub_vertices::Vector{Int64}, sub_flags::Vector{Vector{Int64}}) -> AbstractGKM_subgraph

Return the GKM subgraph induced by the given vertices and flags.
The subgraph includes vertex `sub_vertices[i]` from the supergraph as vertex `i` in the subgraph.
For each vertex `i` in the subgraph, `sub_flags[i]` specifies which flags from the supergraph
vertex `sub_vertices[i]` are included in the subgraph.

The edges of the subgraph are determined by the flags: an edge is included if and only if
both of its corresponding flags are included in the flag subset.

!!! note
    1. This does not check if the result is a valid GKM graph (use `isvalid` for that).
    2. If possible, the subgraph will be endowed with the connection induced from the supergraph.
    3. To create a compact GKM subgraph, it is easier to use `gkm_subgraph_from_vertices` or `gkm_subgraph_from_edges` instead.

# Example
```jldoctest subgr_from_flags
julia> P2 = projective_space(GKM_graph, 2);

julia> # Add standalone flags to make it non-compact
       add_standalone_flag!(P2, 1, gens(P2.M)[1]);

julia> add_standalone_flag!(P2, 2, gens(P2.M)[2]);

julia> add_standalone_flag!(P2, 3, gens(P2.M)[3]);

julia> # Create a subgraph with vertices [1,2] and specific flags
       # Vertex 1 has flags [1, 3] (first edge flag and the standalone flag)
       # Vertex 2 has flags [1, 3] (first edge flag and the standalone flag)
       S = gkm_subgraph_from_flags(P2, [1, 2], [[1, 3], [1, 3]])
GKM subgraph of:
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Standalone flags:
1.3 => (1, 0, 0)
2.3 => (0, 1, 0)
3.3 => (0, 0, 1)
Subgraph:
GKM graph with 2 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
Standalone flags:
1.2 => (1, 0, 0)
2.2 => (0, 1, 0)
```
"""
function gkm_subgraph_from_flags(gkm::AbstractGKM_graph, sub_vertices::Vector{Int64}, sub_flags::Vector{Vector{Int64}}) :: AbstractGKM_subgraph

  @req all(v -> v > 0 && v <= n_vertices(gkm.g), sub_vertices) "Vertex index out of range"
  @req length(sub_vertices) == length(sub_flags) "sub_vertices and sub_flags must have the same length"
  @req length(unique(sub_vertices)) == length(sub_vertices) "sub_vertices has duplicate element"
  if !isvalid(gkm)
    @warn "Creating GKM subgraph of invalid gkm graph. This may result in undefined behavior."
  end

  subnv = length(sub_vertices)

  # Check that flag indices are valid
  if subnv > 0
    for i in 1:subnv
      @req all(f -> f > 0 && f <= valency(gkm), sub_flags[i]) "Flag index out of range for vertex $(sub_vertices[i])"
      @req length(unique(sub_flags[i])) == length(sub_flags[i]) "Duplicate flag indices for vertex $i"
    end
  end

  # Create the underlying graph and determine which edges to include
  g_sub = Graph{Undirected}(subnv)
  # edges_to_add = Tuple{Int64, Int64, AbstractAlgebra.Generic.FreeModuleElem{gkm.weightType}}[]

  for e_super in edges(gkm.g)
    v_super = src(e_super)
    w_super = dst(e_super)

    # Check if both vertices are in the subgraph
    v_sub_idx = findfirst(==(v_super), sub_vertices)
    w_sub_idx = findfirst(==(w_super), sub_vertices)

    if !isnothing(v_sub_idx) && !isnothing(w_sub_idx)
      # Check if the edge flags are both in the subgraph
      flag_at_v = gkm.edge_to_flag_index[e_super]
      flag_at_w = gkm.edge_to_flag_index[reverse(e_super)]

      if flag_at_v in sub_flags[v_sub_idx] && flag_at_w in sub_flags[w_sub_idx]
        # This edge should be in the subgraph
        add_edge!(g_sub, v_sub_idx, w_sub_idx)
        # push!(edges_to_add, (v_sub_idx, w_sub_idx, gkm.w[e_super]))
      end
    end
  end

  # Build the full GKM graph structure manually
  labels = [gkm.labels[sub_vertices[i]] for i in 1:subnv]
  weights_at_vertex = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{gkm.weightType}}}(undef, subnv)
  flag_to_edge = Vector{Vector{Union{Nothing, Edge}}}(undef, subnv)
  edge_to_flag_index = Dict{Edge, Int64}()
  w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{gkm.weightType}}()

  # Build flag structures
  for i in 1:subnv
    v_super = sub_vertices[i]
    num_flags = length(sub_flags[i])
    weights_at_vertex[i] = [gkm.weights_at_vertex[v_super][sub_flags[i][j]] for j in 1:num_flags]
    flag_to_edge[i] = Vector{Union{Nothing, Edge}}(undef, num_flags)

    # Map flags to edges
    for j in 1:num_flags
      flag_idx_super = sub_flags[i][j]
      edge_super = gkm.flag_to_edge[v_super][flag_idx_super]

      if !isnothing(edge_super)
        # This flag corresponds to an edge in the supergraph
        # Check if this edge is in the subgraph
        # v_other_super = src(edge_super) == v_super ? dst(edge_super) : src(edge_super)
        @req src(edge_super) == v_super "Edge-flag relations broken in supergraph. Check with isvalid."
        v_other_super = dst(edge_super)
        v_other_sub_idx = findfirst(==(v_other_super), sub_vertices)

        if !isnothing(v_other_sub_idx) && has_edge(g_sub, i, v_other_sub_idx)
          # This edge is in the subgraph
          e_sub = Edge(i, v_other_sub_idx)
          flag_to_edge[i][j] = e_sub
          edge_to_flag_index[e_sub] = j
          w[e_sub] = gkm.w[edge_super]
        else
          # Edge not in subgraph - this is a standalone flag
          flag_to_edge[i][j] = nothing
        end
      else
        # Standalone flag in supergraph - remains standalone in subgraph
        flag_to_edge[i][j] = nothing
      end
    end
  end

  # Create the subGKM graph object directly
  subGKM = AbstractGKM_graph(
    g_sub,
    labels,
    gkm.M,
    weights_at_vertex,
    edge_to_flag_index,
    flag_to_edge,
    w,
    nothing,  # equivariantCohomology
    nothing,  # curveClasses
    nothing,  # connection
    Dict{CurveClass_type, Array{Any, 3}}(),  # QH_structure_consts
    false  # know_all_QH_structure_consts
  )

  # Initialize the cohomology ring (must be done before initialize!)
  subGKM.equivariantCohomology = _equivariant_cohomology_ring(subGKM)

  res = AbstractGKM_subgraph(gkm, subGKM, sub_vertices, sub_flags)
  _infer_GKM_connection!(res)
  return res
end

@doc"""
    gkm_subgraph_from_edges(gkm::AbstractGKM_graph, edges::Vector{Edge}; include_standalone_flags::Bool=false) -> AbstractGKM_subgraph

Return the GKM subgraph induced by the given edges.

# Arguments
- `gkm`: The supergraph
- `edges`: Vector of edges to include in the subgraph
- `include_all_flags`: If `true`, includes all flags at the subgraph vertices. If `false` (default), only includes edge flags.

!!! note
    1. This does not check if the result is a valid GKM graph (use `isvalid` for that).
    2. If possible, the subgraph will be endowed with the connection induced from the supergraph.
    3. If `include_all_flags=true` the result will have the same valency as the given GKM graph.

# Example
```jldoctest subgr_from_edges
julia> G = projective_space(GKM_graph, 3);

julia> S = gkm_subgraph_from_edges(G, [Edge(1, 2), Edge(2, 3)])
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 3 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 2 => (0, -1, 1, 0)
```
"""
function gkm_subgraph_from_edges(gkm::AbstractGKM_graph, edges::Vector{Edge}; include_all_flags::Bool=false) :: AbstractGKM_subgraph

  # Collect all vertices involved in the given edges
  vDict = Int64[]
  for e in edges
    @req has_edge(gkm.g, e) "Edge $e not found in GKM graph"
    if !(src(e) in vDict)
      push!(vDict, src(e))
    end
    if !(dst(e) in vDict)
      push!(vDict, dst(e))
    end
  end
  sort!(vDict)

  # Build sub_flags: for each vertex, include flags that correspond to the given edges
  # (and optionally all flags if include_all_flags=true)
  subnv = length(vDict)
  sub_flags = Vector{Vector{Int64}}(undef, subnv)

  for i in 1:subnv
    v_super = vDict[i]
    sub_flags[i] = Int64[]

    if include_all_flags
      # Include all flags at this vertex
      append!(sub_flags[i], 1:valency(gkm))
    else
      # Only include flags corresponding to the given edges
      for e in edges
        if src(e) == v_super
          flag_idx = gkm.edge_to_flag_index[e]
          if !(flag_idx in sub_flags[i])
            push!(sub_flags[i], flag_idx)
          end
        elseif dst(e) == v_super
          flag_idx = gkm.edge_to_flag_index[reverse(e)]
          if !(flag_idx in sub_flags[i])
            push!(sub_flags[i], flag_idx)
          end
        end
      end
      sort!(sub_flags[i])
    end
  end

  return gkm_subgraph_from_flags(gkm, vDict, sub_flags)
end

# Return true if the gkm subgraph contains the given edge of the supergraph.
function has_edge(gkmSub::AbstractGKM_subgraph, e::Edge)::Bool
  if !(src(e) in gkmSub.vDict) || !(dst(e) in gkmSub.vDict)
    return false
  end
  sd = indexin([src(e), dst(e)], gkmSub.vDict)
  return has_edge(gkmSub.self.g, Edge(sd[1], sd[2]))
end

function _vertex_preimage(gkmSub::AbstractGKM_subgraph, v::Int64)::Int64
  return indexin([v], gkmSub.vDict)[1]
end

# Return true if the gkm subgraph contains the given vertex of the supergraph.
function has_vertex(gkmSub::AbstractGKM_subgraph, v::Int64)::Bool
  return v in gkmSub.vDict
end

# Return true if the gkm subgraph contains the given vertex label of the supergraph.
function has_vertex(gkmSub::AbstractGKM_subgraph, vertexLabel::String)::Bool
  @req vertexLabel in gkmSub.super.labels "Vertex with label $vertexLabel does not exist."
  v::Int64 = indexin([vertexLabel], gkmSub.super.labels)[1]
  return has_vertex(gkmSub, v)
end

function vertexToSupgraph(gkmSub::AbstractGKM_subgraph, v::Int64)::Int64
  @req v in 1:n_vertices(gkmSub.self.g) "Vertex $v not in subgraph."
  return gkmSub.vDict[v]
end

function edgeToSupergraph(gkmSub::AbstractGKM_subgraph, e::Edge)::Edge
  @req has_edge(gkmSub.self.g, e) "Edge $e not contained in subgraph."
  imE = Edge(vertexToSupgraph(gkmSub, src(e)), vertexToSupgraph(gkmSub, dst(e)))
  return imE
end

@doc raw"""
    is_compatible_with_connection(gkmSub::AbstractGKM_subgraph, con::GKM_connection; printDiagnostics::Bool=true)::Bool

Return `true` if the connection map sends flags contained in the subgraph to flags within the subgraph.
This is necessary for the subgraph to represent a $T$-invariant subspace.
"""
function is_compatible_with_connection(gkmSub::AbstractGKM_subgraph, con::GKM_connection; printDiagnostics::Bool=true)::Bool
  # Check if the connection's valency matches the supergraph's valency
  # This can fail if the supergraph had standalone flags added after the connection was created
  for e_check in edges(gkmSub.super.g)
    if haskey(con.con, e_check) && length(con.con[e_check]) != valency(gkmSub.super)
      printDiagnostics && println("Connection valency mismatch with supergraph (connection may be outdated)")
      return false
    end
  end

  edge_set = collect(edges(gkmSub.self.g))
  for e in Iterators.flatten((edge_set, reverse.(edge_set)))
    eSup = edgeToSupergraph(gkmSub, e)

    # Iterate over actual number of flags at src(e), not uniform valency
    for i in 1:length(gkmSub.flagDict[src(e)])
      # Get the flag index in the supergraph using flagDict
      iSup = gkmSub.flagDict[src(e)][i]

      # Get the connected flag in the supergraph
      jSup = con.con[eSup][iSup]

      # Check if jSup is in flagDict for the destination vertex
      if !any(gkmSub.flagDict[dst(e)] .== jSup)
        printDiagnostics && println("Connection sends flag $i at edge $e to flag $jSup outside subgraph.")
        return false
      end
    end
  end
  return true
end

function Base.show(io::IO, G::AbstractGKM_subgraph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM subgraph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM subgraph with $(n_vertices(G.self.g)) nodes and valency $(valency(G.self))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::AbstractGKM_subgraph)

  println(io, "GKM subgraph of:")
  show(io, MIME"text/plain"(), G.super)
  println(io, "\nSubgraph:")
  show(io, MIME"text/plain"(), G.self)
end

@doc raw"""
    isvalid(gkmsub::AbstractGKM_subgraph; printDiagnostics::Bool = true) -> Bool

Return true if the given GKM subgraph is valid. This holds if and only if all of the following hold:
  1. The supergraph and subgraph are both valid GKM graphs of the same character group
  2. The subgraph is mathematically a subgraph of the supergraph
  3. The edge weights of the subgraph match that of the supergraph
  4. The vertex labels of the subgraph and the supergraph match
  5. The flag weights of the subgraph match the corresponding flags in the supergraph
  6. flagDict correctly maps flags from subgraph to supergraph
!!! warning
    If a connection for the supergraph is set, this does not check if it is compatible with the subgraph.
    Use `is_compatible_with_connection()` for this.
"""
function isvalid(gkmsub::AbstractGKM_subgraph; printDiagnostics::Bool = true)::Bool
  if !isvalid(gkmsub.super; printDiagnostics)
    printDiagnostics && println("GKM-Supergraph is invalid")
    return false
  elseif !isvalid(gkmsub.self; printDiagnostics)
    printDiagnostics && println("Sub-GKM-graph is invalid as GKM graph")
    return false
  elseif gkmsub.self.M != gkmsub.super.M
    printDiagnostics && println("GKM parent and subgraph don't have the same character group")
    return false
  end

  parentVertices = 1:n_vertices(gkmsub.super.g)
  for v in gkmsub.vDict
    if !(v in parentVertices)
      printDiagnostics && println("Vertex $v not in parent GKM graph")
      return false
    end
  end
  for e in edges(gkmsub.self.g)
    targetEdge = Edge(gkmsub.vDict[src(e)], gkmsub.vDict[dst(e)])
    if !has_edge(gkmsub.super.g, targetEdge)
      printDiagnostics && println("Edge $e gets mapped to non-existent edge $targetEdge in parent GKM graph")
      return false
    elseif gkmsub.self.w[e] != gkmsub.super.w[targetEdge]
      printDiagnostics && println("Weights of $e and its image $targetEdge in the parent GKM graph don't match")
      return false
    end
  end

  for v in 1:n_vertices(gkmsub.self.g)
    if gkmsub.self.labels[v] != gkmsub.super.labels[gkmsub.vDict[v]]
      printDiagnostics && println("Label of vertex $v disagrees in subgraph and supergraph.")
      return false
    end
  end

  # Check flag consistency using flagDict
  if length(gkmsub.flagDict) != n_vertices(gkmsub.self.g)
    printDiagnostics && println("flagDict length doesn't match number of vertices in subgraph")
    return false
  end

  for v_sub in 1:n_vertices(gkmsub.self.g)
    v_super = gkmsub.vDict[v_sub]

    # Check that all flag indices are valid
    for flag_idx_sub in 1:valency(gkmsub.self)
      flag_idx_super = gkmsub.flagDict[v_sub][flag_idx_sub]

      if flag_idx_super < 1 || flag_idx_super > valency(gkmsub.super)
        printDiagnostics && println("Invalid flag index in flagDict at vertex $v_sub")
        return false
      end

      # Check that the flag weights match
      weight_sub = gkmsub.self.weights_at_vertex[v_sub][flag_idx_sub]
      weight_super = gkmsub.super.weights_at_vertex[v_super][flag_idx_super]

      if weight_sub != weight_super
        printDiagnostics && println("Flag weight mismatch at vertex $v_sub, flag $flag_idx_sub")
        return false
      end
    end

    # Check for duplicate flag indices
    if length(unique(gkmsub.flagDict[v_sub])) != length(gkmsub.flagDict[v_sub])
      printDiagnostics && println("Duplicate flag indices in flagDict at vertex $v_sub")
      return false
    end
  end

  return true
end