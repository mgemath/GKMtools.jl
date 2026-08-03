# This file is part of GKMtools.jl, licensed under the MIT License (MIT).

# R is the type of the character lattice (typically ZZ or QQ module)
# V is the type of the vertex labels (typically String)
# F is the type of the flag weights (typically ToricFlagWeight or OrbifoldToricFlagWeight)

struct GKMCombinatorialData{R, V, F}
  
  ###########################################################################
  # CORE COMBINATORIAL DATA
  ###########################################################################

  g::Graph{Undirected}

  # Character lattice (typically ZZ or QQ module)
  M::AbstractAlgebra.Generic.FreeModule{R}

  # Vertex labels
  labels::Vector{V}

  # flags[v] = list of flags weights at vertex v
  flags::Vector{Vector{F}}

  # edge -> (flag index at source, flag index at target)
  edge_flags::Dict{Edge, Tuple{Int,Int}}
end

# --- GKMCombinatorialData ---

function Base.show(io::IO, data::GKMCombinatorialData)
  print(io, "GKM Combinatorial Data with $(nv(data.g)) vertices and $(ne(data.g)) edges")
end

function Base.show(io::IO, ::MIME"text/plain", data::GKMCombinatorialData)
  println(io, "GKM Combinatorial Data")
  println(io, "  Graph: ", data.g)
  println(io, "  Character Lattice: ", data.M)
  println(io, "  Vertex Labels: ", get_string.(data.labels))
  print(io, "  Flag weights defined at $(length(data.flags)) vertices")
end
