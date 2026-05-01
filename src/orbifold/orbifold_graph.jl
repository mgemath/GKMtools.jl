struct OrbifoldGKMGraph{R} <: AbstractGKMGraph{R}
  base::GKMGraph{R}

  vertex_isotropy::Vector{Vector{Int}}
  edge_multiplicity::Dict{Edge,Int}
end

# Interface delegation

graph(g::OrbifoldGKMGraph) = graph(g.base)
lattice(g::OrbifoldGKMGraph) = lattice(g.base)
flags(g::OrbifoldGKMGraph, v::Int) = flags(g.base, v)

function weight(g::OrbifoldGKMGraph, e::Edge)
  return g.edge_multiplicity[e] * weight(g.base, e)
end

num_vertices(g::OrbifoldGKMGraph) = num_vertices(g.base)
num_edges(g::OrbifoldGKMGraph) = num_edges(g.base)