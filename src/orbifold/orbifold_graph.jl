# This file is part of GKMtools.jl, licensed under the MIT License (MIT).

# Isotropy data for vertices and edge multiplicities for orbifold GKM graphs
# For each vertex, we store the invariants of the isotropy group and the weights of the representation on the tangent space. 
# For each edge, we store the multiplicity (the order of the isotropy group along that edge).  
struct VertexIsotropyData
    invariants::Vector{Int}   # group structure of the isotropy group at the vertex
    weights::Matrix{Int}      # full tangent representation
end

struct FlagIsotropyData
    invariants::Vector{Int}   # group structure of the isotropy group at the edge
    weights::Matrix{Int}      # full tangent representation
end

struct OrbifoldGKMGraph{R} <: AbstractGKMGraph{R}
  base::GKMGraph{R}
  vertex_isotropy::Vector{VertexIsotropyData}
  flag_isotropy::Dict{Flag, FlagIsotropyData}
end

# Interface delegation

graph(g::OrbifoldGKMGraph) = graph(g.base)
lattice(g::OrbifoldGKMGraph) = lattice(g.base)
flags(g::OrbifoldGKMGraph, v::Int) = flags(g.base, v)

# function weight(g::OrbifoldGKMGraph, e::Edge)
#   return g.edge_multiplicity[e] * weight(g.base, e)
# end


num_vertices(g::OrbifoldGKMGraph) = num_vertices(g.base)
num_edges(g::OrbifoldGKMGraph) = num_edges(g.base)

function Base.show(io::IO, ogkm::OrbifoldGKMGraph{R}) where R
    # Multi-line format is often better for complex nested structs
    println(io, "OrbifoldGKMGraph{$R}:")
    print(io, "  Base Graph: ")
    show(io, ogkm.base) 
    println(io)
    print(io, "  Isotropy: $(length(ogkm.vertex_isotropy)) vertices, ")
    print(io, "$(length(ogkm.flag_isotropy)) flags")
end