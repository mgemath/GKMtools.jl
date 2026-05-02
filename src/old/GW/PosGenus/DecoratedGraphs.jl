# Create a decorated graph corresponding to the given GKM graph from the given data.
# If the graph has a single vertex and total genus zero, it must have at least three marked points.
# If it has a single vertex and total genus one, it must have at least one marked point.
# This is a generalization of decoratedTree(...) to postive genus.
function decoratedGraph(
  gkm::AbstractGKM_graph,
  g::Graph,
  vDict::Union{Vector{Int}, Tuple{Vararg{Int}}},
  edgeMult::Dict{Edge, Vector{Int}},
  marks::Vector{Int},
  genus::Vector{Int};
  check::Bool=true)::GW_decorated_graph

  nv = n_vertices(g)

  if check && nv == 1 && degree(g, 1) == 0 && genus[1]==0
    @req length(marks) >= 3 "Single vertex decorated graph of genus zero must have at least three marked points"
  elseif check && nv == 1 && degree(g, 1) == 0 && genus[1]==1
    @req length(marks) >= 1 "Single vertex decorated graph of genus one must have at least one marked point"
  end

  for e in collect(keys(edgeMult))
    edgeMult[reverse(e)] = edgeMult[e]
  end

  if check
    @req nv == length(vDict) "g and vDict have different lengths"
    @req all([has_edge(gkm.g, Edge(vDict[src(e)], vDict[dst(e)])) for e in edges(g)]) "image of edge does not exist in GKM graph"
    @req all(e -> all(m -> 0 < m, edgeMult[e]), edges(g)) "non-positive edge multiplicity"
    @req all(g -> g >= 0, genus) "Genus markings of vertices must be non-negative"
  end

  return GW_decorated_graph(gkm, g, vDict, edgeMult, marks, genus)
end

# Return the image of the edge e of the decorated graph g in the corresponding GKM graph.
function imageOf(e::Edge, t::GW_decorated_graph)::Edge
  return Edge(t.vDict[src(e)], t.vDict[dst(e)])
end

# Return the image of the vertex v of the decorated graph t in the corresponding GKM graph
function imageOf(v::Int, t::GW_decorated_graph)::Int
  return t.vDict[v]
end

function edgeMult(e::Edge, dt::GW_decorated_graph)::Vector{Int}
  if e in keys(dt.edgeMult)
    return dt.edgeMult[e]
  elseif reverse(e) in keys(dt.edgeMult)
    return dt.edgeMult[reverse(e)]
  end
  @req false "edge has no multiplicity assigned in decorated graph"
end

function valency(v::Int, dg::GW_decorated_graph)::Int
  return sum(n -> length(dg.edgeMult[Edge(v, n)]), all_neighbors(dg.g, v))
end