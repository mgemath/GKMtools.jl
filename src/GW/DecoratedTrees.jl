# Create a decorated tree corresponding to the given GKM graph from the given data.
# If the tree has a single vertex, it must have at least three marked points.
# Warning: This does not check if the given tree is actually a tree.
function decoratedTree(
  gkm::AbstractGKM_graph,
  tree::Graph,
  vDict::Union{Vector{Int}, Tuple{Vararg{Int}}},
  edgeMult::Dict{Edge, Int},
  marks::Vector{Int};
  check::Bool=true)::GW_decorated_tree

  nv = n_vertices(tree)

  if check && nv == 1 && degree(tree, 1) == 0
    @req length(marks) >= 3 "Single vertex tree must have at least three marked points"
  end

  for e in collect(keys(edgeMult))
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
function imageOf(v::Int, t::GW_decorated_tree)::Int
  return t.vDict[v]
end

function edgeMult(e::Edge, dt::GW_decorated_tree)::Int
  if e in keys(dt.edgeMult)
    return dt.edgeMult[e]
  elseif reverse(e) in keys(dt.edgeMult)
    return dt.edgeMult[reverse(e)]
  end
  @req false "edge has no multiplicity assigned in decorated tree"
end

# This function implements Theorem 4.3 from Giosue's "Enumeration of rational contact curves via torus actions".
# It assumes that dt is a tree with map to some projective space.
# TODO: Giosue, please check if t1, t2, ... vs -t1, -t2, -t3 is correct here.
function _incidency(dt::GW_decorated_tree, r::Int64)::QQMPolyRingElem
  R = dt.gkm.equivariantCohomology.coeffRing
  res = zero(R)
  l = gens(R)
  for e in edges(dt.tree)
    eRes = zero(R)
    for t in 0:r
      eRes += l[src(e)]^(t) * l[dst(e)]^(r-t)
    end
    res += eRes * dt.edgeMult[e]
  end
  return res
end