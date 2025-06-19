@doc"""
    gkm_subgraph_from_vertices(gkm::GKM_graph, vertices::Vector{Int64}) -> GKM_subgraph

Return the GKM subgraph induced by the given vertices.
!!! note
    1. This does not check if the result is a valid GKM graph (use may use `isvalid` for that).
    2. If possible, the subgraph will be endowed with the connection induced from the supergraph.

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

julia> S.Domain
GKM graph with 2 nodes, valency 1 and axial function:
3 -> 2 => (0, -1, 1, 0)

julia> S.Codomain
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

```
"""
function gkm_subgraph_from_vertices(gkm::GKM_graph, vertices::Vector{Int64})::GKM_subgraph

  @req all(v -> v>0, vertices) "Vertex index must be positive"

  return _gkm_subgraph_from_vertices(gkm, unique(sort(vertices)))
end

function _gkm_subgraph_from_vertices(gkm::GKM_graph, vDict::Vector{Int64})::GKM_subgraph

  subnv = length(vDict)
  labels = [gkm.labels[vDict[i]] for i in 1:subnv]
  subGKM = gkm_graph(Graph{Undirected}(subnv), labels, gkm.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{typeof(_get_weight_type(gkm))}}(); checkLabels=false) # use the same character lattice
  
  for e in edges(gkm.g)
    if src(e) in vDict && dst(e) in vDict
      sd = indexin([src(e), dst(e)], vDict)
      add_edge!(subGKM, sd[1], sd[2], gkm.w[e])
    end
  end
  res = GKM_subgraph(subGKM, gkm, vDict)
  _infer_GKM_connection!(res)
  return res
end


# If the supergraph's connection is compatible with the subgraph, infer it to the subgraph and return true.
# Else return false. False is also returned if the subgraph's connection is already set.
function _infer_GKM_connection!(gkmSub::GKM_subgraph)::Bool

  con = get_connection(gkmSub.Codomain)
  if !isnothing(con) && is_compatible_with_connection(gkmSub, con; printDiagnostics=false)

    oldCon = con.con
    oldA = con.a
    newCon = Dict{Tuple{Edge, Edge}, Edge}()
    newA = Dict{Tuple{Edge, Edge}, ZZRingElem}()
    vDict = gkmSub.vDict

    for v in 1:n_vertices(gkmSub.Domain.g)
      for w in all_neighbors(gkmSub.Domain.g, v)
        for u in all_neighbors(gkmSub.Domain.g, v)
          e = Edge(v, w)
          ei = Edge(v, u)
          eSup = Edge(vDict[v], vDict[w])
          eiSup = Edge(vDict[v], vDict[u])
          epiSup = oldCon[(eSup, eiSup)]
          epi = Edge(w, _vertex_preimage(gkmSub, dst(epiSup)))
          newCon[(e, ei)] = epi
          newA[(e, ei)] = oldA[(eSup, eiSup)]
        end
      end
    end

    newConObj = GKM_connection(newCon, newA)
    set_connection!(gkmSub.Domain, newConObj)
    return true
  end
  return false
end

@doc"""
    gkm_subgraph_from_vertices(gkm::GKM_graph, vertexLabels::Vector{String}) -> GKM_subgraph

As before, but the vertices are given by their labels.
"""
function gkm_subgraph_from_vertices(gkm::GKM_graph, vertexLabels::Vector{String})::GKM_subgraph

  @req all(l -> l in gkm.labels, vertexLabels) "Label not found"

  vertices::Vector{Int64} = indexin(vertexLabels, gkm.labels) # need to specify Vector{Int64} as indexin returns vector of Union{Nothing, Int64}.
  res = gkm_subgraph_from_vertices(gkm, vertices)
  _infer_GKM_connection!(res)
  return res
end

@doc"""
    gkm_subgraph_from_edges(gkm::GKM_graph, edges::Vector{Edge}) -> GKM_subgraph
    
Return the GKM subgraph induced by the given edges.
!!! note
    1. This does not check if the result is a valid GKM graph (use `isvalid` for that).
    2. If possible, the subgraph will be endowed with the connection induced from the supergraph.

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
function gkm_subgraph_from_edges(gkm::GKM_graph, edges::Vector{Edge})::GKM_subgraph
  
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

  subnv = length(vDict)
  labels = [gkm.labels[vDict[i]] for i in 1:subnv]
  subGKM = gkm_graph(Graph{Undirected}(subnv), labels, gkm.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{gkm.weightType}}(); checkLabels=false) # use the same character lattice
  
  for e in edges
    sd = indexin([src(e), dst(e)], vDict)
    add_edge!(subGKM, sd[1], sd[2], gkm.w[e])
  end
  res = GKM_subgraph(subGKM, gkm, vDict)
  _infer_GKM_connection!(res)
  return res
end

# Return true if the gkm subgraph contains the given edge of the supergraph.
function has_edge(gkmSub::GKM_subgraph, e::Edge)::Bool
  if !(src(e) in gkmSub.vDict) || !(dst(e) in gkmSub.vDict)
    return false
  end
  sd = indexin([src(e), dst(e)], gkmSub.vDict)
  return has_edge(gkmSub.Domain.g, Edge(sd[1], sd[2]))
end

function _vertex_preimage(gkmSub::GKM_subgraph, v::Int64)::Int64
  return indexin([v], gkmSub.vDict)[1]
end

# Return true if the gkm subgraph contains the given vertex of the supergraph.
function has_vertex(gkmSub::GKM_subgraph, v::Int64)::Bool
  return v in gkmSub.vDict
end

# Return true if the gkm subgraph contains the given vertex label of the supergraph.
function has_vertex(gkmSub::GKM_subgraph, vertexLabel::String)::Bool
  @req vertexLabel in gkmSub.Codomain.labels "Vertex with label $vertexLabel does not exist."
  v::Int64 = indexin([vertexLabel], gkmSub.Codomain.labels)[1]
  return has_vertex(gkmSub, v)
end

function vertexToSupgraph(gkmSub::GKM_subgraph, v::Int64)::Int64
  @req v in 1:n_vertices(gkmSub.Domain.g) "Vertex $v not in subgraph."
  return gkmSub.vDict[v]
end

function edgeToSupergraph(gkmSub::GKM_subgraph, e::Edge)::Edge
  @req has_edge(gkmSub.Domain.g, e) "Edge $e not contained in subgraph."
  imE = Edge(vertexToSupgraph(gkmSub, src(e)), vertexToSupgraph(gkmSub, dst(e)))
  return imE
end

@doc raw"""
    is_compatible_with_connection(gkmSub::GKM_subgraph, con::GKM_connection; printDiagnostics::Bool=true)::Bool

Return `true` if the connection map sends edge pairs contained in the subgraph to an edge of the subgraph.
This is necessary for the subgraph to represent a $T$-invariant subspace.
"""
function is_compatible_with_connection(gkmSub::GKM_subgraph, con::GKM_connection; printDiagnostics::Bool=true)::Bool
  nvsub = n_vertices(gkmSub.Domain.g)
  for v in 1:nvsub
    for w in 1:nvsub
      e = Edge(v, w)
      if !has_edge(gkmSub.Domain.g, e)
        continue
      end
      for u in all_neighbors(gkmSub.Domain.g, v)
        ei = Edge(v, u)
        epi = con.con[(edgeToSupergraph(gkmSub, e), edgeToSupergraph(gkmSub, ei))]
        if !has_edge(gkmSub, epi)
          printDiagnostics && println("Connection sends image of ($e, $ei) in supergraph to outside the subgraph.")
          return false
        end
      end
    end
  end
  return true
end

function Base.show(io::IO, G::GKM_subgraph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM subgraph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM subgraph with $(n_vertices(G.Domain.g)) nodes and valency $(valency(G.Domain))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::GKM_subgraph)

  println(io, "GKM subgraph of:")
  show(io, MIME"text/plain"(), G.Codomain)
  println(io, "\nSubgraph:")
  show(io, MIME"text/plain"(), G.Domain)
end

@doc raw"""
    isvalid(gkmsub::GKM_subgraph; printDiagnostics::Bool = true) -> Bool

Return true if the given GKM subgraph is valid. This holds if and only if all of the following hold:
  1. The supergraph and subgraph are both valid GKM GKMsubgraphs of the same character group
  2. The subgraph is mathematically a subgraph of the supergraph
  3. The edge weights of the subgraph match that of the supergraph
  4. The vertex labels of the subgraph and the supergraph match.
!!! warning
    If a connection for the supergraph is set, this does not check if it is compatible with the subgraph.
    Use `is_compatible_with_connection()` for this.
"""
function isvalid(gkmsub::GKM_subgraph; printDiagnostics::Bool = true)::Bool
  if !isvalid(gkmsub.Codomain; printDiagnostics)
    printDiagnostics && println("GKM-Supergraph is invalid")
    return false
  elseif !isvalid(gkmsub.Domain; printDiagnostics)
    printDiagnostics && println("Sub-GKM-graph is invalid as GKM graph")
    return false
  elseif gkmsub.Domain.M != gkmsub.Codomain.M
    printDiagnostics && println("GKM parent and subgraph don't have the same character group")
    return false
  end
  
  parentVertices = 1:n_vertices(gkmsub.Codomain.g)
  for v in gkmsub.vDict
    if !(v in parentVertices)
      printDiagnostics && println("Vertex $v not in parent GKM graph")
      return false
    end
  end
  for e in edges(gkmsub.Domain.g)
    targetEdge = Edge(gkmsub.vDict[src(e)], gkmsub.vDict[dst(e)])
    if !has_edge(gkmsub.Codomain.g, targetEdge)
      printDiagnostics && println{"Edge $e in gets mapped to non-existent edge $targetEdge in parent GKM graph"}
      return false
    elseif gkmsub.Domain.w[e] != gkmsub.Codomain.w[targetEdge]
      printDiagnostics && println("Weights of $e and its image $targetEdge in the parent GKM graph don't match")
      return false
    end
  end
  
  for v in 1:n_vertices(gkmsub.Domain.g)
    if gkmsub.Domain.labels[v] != gkmsub.Codomain.labels[gkmsub.vDict[v]]
      printDiagnostics && println("Label of vertex $v disagrees in subgraph and supergraph.")
      return false
    end
  end
  return true
end