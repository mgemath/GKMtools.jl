struct GW_decorated_graph{G<:AbstractGKMGraph,C}
  gkm::G
  g::Graph{Undirected}
  vDict::Vector{Int}
  edgeMult::Dict{Edge,Vector{Int}}
  marks::Vector{Int}
  genus::Vector{Int}
  context::C
end

function decoratedGraph(gkm::AbstractGKMGraph, g::Graph{Undirected}, vDict,
  edgeMult::Dict{Edge,Vector{Int}}, marks::Vector{Int}, genus::Vector{Int},
  context::GWClassEvaluationContext; check::Bool=true)
  if check
    @req nv(g) == length(vDict) "g and vDict have different lengths"
    @req all(e -> has_edge(graph(gkm), Edge(vDict[src(e)], vDict[dst(e)])), edges(g)) "image edge does not exist"
    @req all(e -> all(>(0), edgeMult[e]), edges(g)) "non-positive edge multiplicity"
    @req all(>=(0), genus) "vertex genera must be non-negative"
  end
  return GW_decorated_graph(gkm, g, collect(vDict), edgeMult, marks, genus, context)
end

imageOf(e::Edge, graph::GW_decorated_graph) = Edge(graph.vDict[src(e)], graph.vDict[dst(e)])
imageOf(v::Int, graph::GW_decorated_graph) = graph.vDict[v]

function edgeMult(e::Edge, graph::GW_decorated_graph)
  haskey(graph.edgeMult, e) && return graph.edgeMult[e]
  haskey(graph.edgeMult, reverse(e)) && return graph.edgeMult[reverse(e)]
  @req false "edge has no multiplicity assigned in decorated graph"
end

valency(v::Int, graph::GW_decorated_graph) =
  sum(n -> length(edgeMult(Edge(v, n), graph)), all_neighbors(graph.g, v))
