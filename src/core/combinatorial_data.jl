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
  println(io, "  Flag weights defined at $(length(data.flags)) vertices")
end

graph(G::AbstractGKMGraph) = G.core.g
num_vertices(G::AbstractGKMGraph) = nv(G.core.g)
num_edges(G::AbstractGKMGraph) = ne(G.core.g)
lattice(G::AbstractGKMGraph) = G.core.M
labels(G::AbstractGKMGraph) = G.core.labels
flags(G::AbstractGKMGraph, v::Int) = G.core.flags[v]

edges(G::AbstractGKMGraph) = edges(G.core.g)
vertices(G::AbstractGKMGraph) = vertices(G.core.g)

label(G::AbstractGKMGraph, v::Int) = G.core.labels[v].label
degree(G::AbstractGKMGraph, v::Int) = degree(G.core.g, v)

function valency(G::AbstractGKMGraph)
  return _valency(G, check = false)
end

function _valency(G::AbstractGKMGraph; check::Bool = true)
  if check
    for v in vertices(G)
      if length(flags(G, 1)) != length(flags(G, v))
        error("Valency check failed at vertex $v: number of flags does not match degree")
      end
    end
  end

  return length(flags(G, 1))
end

function is_compact(G::AbstractGKMGraph)
  val = valency(G)
  for v in vertices(G)
    if length(flags(G, v)) != val
      return false
    end
  end
  return true
end

function compact_flags(G::AbstractGKMGraph, v::Int)
  
  if degree(G, v) == length(flags(G, v))
    return collect(1:length(flags(G, v)))
  end

  ans = Vector{Int}(undef, degree(G, v))
  index = 1
  for e in edges(G)
    if src(e) == v
      i, _ = G.edge_flags[e]
      ans[index] = i
      index += 1
    elseif dst(e) == v
      _, j = G.edge_flags[e]
      ans[index] = j
      index += 1
    end
  end
  
  return sort(ans)
end

function weight(G::AbstractGKMGraph{R, V, F}, e::Edge) where {R, V, F}
  _e = e
  sign = 1
  if !haskey(G.core.edge_flags, e)
    # Try the opposite edge for undirected graphs
    _e = Edge(dst(e), src(e))
    sign = -1
    
    if !haskey(G.core.edge_flags, _e)
      error("Edge $e not found in edge_flags")
    end
  end
  i, _ = G.core.edge_flags[_e]
  v = src(_e)
  return sign * G.core.flags[v][i].weight
end