# Create a decorated tree corresponding to the given GKM graph from the given data.
# If the tree has a single vertex, it must have at least three marked points.
# Warning: This does not check if the given tree is actually a tree.
struct GW_decorated_tree
  gkm::AbstractGKM_graph
  tree::Graph
  vDict::Union{Vector{Int64}} # map vertices of tree to vertices of gkm.g
  edgeMult::Dict{Edge, Int} # each edge of tree has a non-negative multiplicity
  marks::Vector{Int} # vector of marked vertices of the tree 
  
  function GW_decorated_tree(
    gkm::AbstractGKM_graph,
    tree::Graph{Undirected},
    vDict::Union{Vector{Int64}},
    edgeMult::Dict{Edge, Int},
    marks::Vector{Int}
  )
    return new(gkm, tree, vDict, edgeMult, marks)
  end
end

function decoratedTree(
  gkm::GKM_graph,
  tree::Graph,
  vDict::Vector{Int64},
  edgeMult::Dict{Edge, Int64},
  marks::Vector{Int64};
  check::Bool=true)::GW_decorated_tree

  nv = n_vertices(tree)

  if check && nv == 1 && degree(tree, 1) == 0
    @req length(marks) >= 3 "Single vertex tree must have at least three marked points"
  end

  for e in keys(edgeMult)
    edgeMult[reverse(e)] = edgeMult[e]
  end

  if check
    @req nv == length(vDict) "tree and vDict have different lengths"
    @req all([has_edge(gkm.g, Edge(vDict[src(e)], vDict[dst(e)])) for e in edges(tree)]) "image of tree edge does not exist in GKM graph"
    @req all([0 < edgeMult[e] for e in edges(tree)]) "non-positive edge multiplicity"
  end

  return GW_decorated_tree(gkm, tree, vDict, edgeMult, marks)
end

# Return the image of the edge e of the decorated tree t in the corresponding GKM graph.
function imageOf(e::Edge, t::GW_decorated_tree)::Edge
  return Edge(t.vDict[src(e)], t.vDict[dst(e)])
end

# Return the image of the vertex v of the decorated tree t in the corresponding GKM graph
function imageOf(v::Int64, t::GW_decorated_tree)::Int64
  return t.vDict[v]
end

function edgeMult(e::Edge, dt::GW_decorated_tree)::Int64
  if e in keys(dt.edgeMult)
    return dt.edgeMult[e]
  elseif reverse(e) in keys(dt.edgeMult)
    return dt.edgeMult[reverse(e)]
  end
  @req false "edge has no multiplicity assigned in decorated tree"
end

