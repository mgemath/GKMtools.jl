"""
    _calculate_graph_cycles_via_trees(G, edgeList, edgeToGenIndex)

Return a basis of the integral first homology of the graph underlying `G`.

A breadth-first spanning tree is constructed in each connected component.
Every edge outside the resulting spanning forest determines one fundamental
cycle. Each cycle is returned as an ordered list of oriented edges: the target
of every edge is the source of the next edge, including cyclically from the
last edge back to the first.
"""
function _calculate_graph_cycles_via_trees(
  G::AbstractGKM_graph,
  edgeList::Vector{Edge},
  edgeToGenIndex::Dict{Edge, Int64}
)::Vector{Vector{Edge}}
  nVertices = n_vertices(G.g)
  nEdges = n_edges(G.g)

  @assert length(edgeList) == nEdges

  visited = falses(nVertices)
  parentVertex = zeros(Int, nVertices)
  depth = zeros(Int, nVertices)
  isTreeEdge = falses(nEdges)

  # A vector with a moving head is enough for all BFS queues and avoids an
  # additional queue dependency.
  queue = Vector{Int}(undef, nVertices)
  nComponents = 0

  for root in 1:nVertices
    visited[root] && continue
    nComponents += 1
    visited[root] = true

    head = 1
    tail = 1
    queue[tail] = root

    while head <= tail
      v = queue[head]
      head += 1

      for w in all_neighbors(G.g, v)
        visited[w] && continue

        eIndex = edgeToGenIndex[Edge(v, w)]
        visited[w] = true
        parentVertex[w] = v
        depth[w] = depth[v] + 1
        isTreeEdge[eIndex] = true

        tail += 1
        queue[tail] = w
      end
    end
  end

  cycles = Vector{Vector{Edge}}()
  sizehint!(cycles, nEdges - nVertices + nComponents)

  for (edgeIndex, e) in enumerate(edgeList)
    isTreeEdge[edgeIndex] && continue

    # The non-tree edge is traversed from src(e) to dst(e), so close the
    # cycle by following tree edges from dst(e) back to src(e). The two
    # vertices are raised to equal depth and then together to their LCA.
    u = src(e)
    v = dst(e)
    cycle = Edge[e]
    uBranch = Edge[]
    sizehint!(cycle, depth[u] + depth[v] + 1)
    sizehint!(uBranch, depth[u])

    while depth[v] > depth[u]
      p = parentVertex[v]
      push!(cycle, Edge(v, p))
      v = p
    end

    while depth[u] > depth[v]
      p = parentVertex[u]
      push!(uBranch, Edge(p, u))
      u = p
    end

    while u != v
      pv = parentVertex[v]
      push!(cycle, Edge(v, pv))
      v = pv

      pu = parentVertex[u]
      push!(uBranch, Edge(pu, u))
      u = pu
    end

    # `uBranch` was found from the non-tree edge's source towards the LCA,
    # whereas the oriented loop traverses it from the LCA back to the source.
    for i in length(uBranch):-1:1
      push!(cycle, uBranch[i])
    end

    push!(cycles, cycle)
  end

  return cycles
end
