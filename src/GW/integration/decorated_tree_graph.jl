struct GW_decorated_tree{G<:AbstractGKMGraph,C}
  gkm::G
  tree::Graph{Undirected}
  vDict::Vector{Int} # map vertices of tree to vertices of gkm.g
  edgeMult::Dict{Edge, Int} # each edge of tree has a non-negative multiplicity
  marks::Vector{Int} # vector of marked vertices of the tree
  context::C
end

function decoratedTree(gkm::AbstractGKMGraph, tree::Graph, vDict,
  edgeMult::Dict{Edge,Int}, marks::Vector{Int}, context::GWClassEvaluationContext; check::Bool=true)
  check && @req nv(tree) == length(vDict) "tree and vDict have different lengths"
  check && @req all(e -> has_edge(graph(gkm), Edge(vDict[src(e)], vDict[dst(e)])), edges(tree)) "image edge does not exist"
  check && @req all(e -> edgeMult[e] > 0, edges(tree)) "non-positive edge multiplicity"
  return GW_decorated_tree(gkm, tree, collect(vDict), edgeMult, marks, context)
end

imageOf(e::Edge, tree::GW_decorated_tree) = Edge(tree.vDict[src(e)], tree.vDict[dst(e)])
imageOf(v::Int, tree::GW_decorated_tree) = tree.vDict[v]

function edgeMult(e::Edge, tree::GW_decorated_tree)
  haskey(tree.edgeMult, e) && return tree.edgeMult[e]
  haskey(tree.edgeMult, reverse(e)) && return tree.edgeMult[reverse(e)]
  @req false "edge has no multiplicity assigned in decorated tree"
end
function Euler_inv(dt::GW_decorated_tree, t::Vector{T}; check_degree::Bool=false) where T
  result = one(t[1])
  for v in vertices(dt.tree)
    val = degree(dt.tree, v)
    euler = _euler_class(dt.gkm, imageOf(v, dt), t)
    result = val >= 1 ? result * euler^(val - 1) : result // euler
    inverse_weights = zero(t[1])
    for neighbor in all_neighbors(dt.tree, v)
      edge = Edge(v, neighbor)
      edge_weight = _weight_class(dt.gkm, imageOf(edge, dt), t) // edgeMult(edge, dt)
      result = result // edge_weight
      inverse_weights += 1 // edge_weight
    end
    exponent = val - 3 + count(==(v), dt.marks)
    result = exponent >= 0 ? result * inverse_weights^exponent : result // inverse_weights^(-exponent)
  end
  return result
end

function _b(u, w, a)
  result = one(u) // one(u)
  if a >= 0
    for j in 0:a
      result = result // (w - j*u)
    end
  else
    for j in 1:(-a-1)
      result *= w + j*u
    end
  end
  return result
end

function _source_flag_index(G::AbstractGKMGraph, e::Edge)
  data = core(G)
  return haskey(data.edge_flags, e) ? first(data.edge_flags[e]) : last(data.edge_flags[reverse(e)])
end

function _h(e::Edge, d::Int, con::AbstractGKMConnection, G::AbstractGKMGraph, t::Vector{T}) where T
  edge_weight = _weight_class(G, e, t)
  result = one(t[1]) * (-1)^d * ZZ(d)^(2d) // (factorial(ZZ(d))^2 * edge_weight^(2d))
  edge_flag = _source_flag_index(G, e)
  for i in 1:valency(G)
    i == edge_flag && continue
    result *= _b(edge_weight // d, _flag_weight_class(G, src(e), i, t), coefficients(con)[e][i] * d)
  end
  return result
end
