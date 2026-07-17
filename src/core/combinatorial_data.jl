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

core(G::AbstractGKMGraph) = G.core
graph(data::GKMCombinatorialData) = data.g
graph(G::AbstractGKMGraph) = graph(core(G))

num_vertices(G::AbstractGKMGraph) = nv(graph(G))
num_edges(G::AbstractGKMGraph) = ne(graph(G))

lattice(data::GKMCombinatorialData) = data.M
lattice(G::AbstractGKMGraph) = lattice(core(G))

labels(data::GKMCombinatorialData) = data.labels
labels(G::AbstractGKMGraph) = labels(core(G))

flags(data::GKMCombinatorialData, v::Int) = data.flags[v]
flags(G::AbstractGKMGraph, v::Int) = flags(core(G), v)

edges(G::AbstractGKMGraph) = edges(graph(G))

vertices(core::GKMCombinatorialData) = vertices(graph(core))
vertices(G::AbstractGKMGraph) = vertices(core(G))

label(G::AbstractGKMGraph, v::Int) = core(G).labels[v].label
degree(G::AbstractGKMGraph, v::Int) = degree(graph(G), v)

rank_torus(core::GKMCombinatorialData) = rank(core.M)
rank_torus(G::AbstractGKMGraph) = rank_torus(core(G))

function find_vertex_index(Vertexlabel::String, G::AbstractGKMGraph)
  index = 0
  for i in 1:num_vertices(G)
    if Vertexlabel == label(G, i)
      index = i
      break
    end
  end
  @assert (index > 0) "label not found"
  return index
end

function valency(G)
  return _valency(G, check = false)
end

function _valency(G; check::Bool = true)
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

function weight(core::GKMCombinatorialData, v::Int, i::Int)
  return core.flags[v][i].weight
end

function weight(G::AbstractGKMGraph, v::Int, i::Int)
  return weight(core(G), v, i)
end

function weight(core::GKMCombinatorialData, e::Edge)
  _e = e
  sign = 1
  if !haskey(core.edge_flags, e)
    # Try the opposite edge for undirected graphs
    _e = Edge(dst(e), src(e))
    sign = -1
    
    if !haskey(core.edge_flags, _e)
      error("Edge $e not found in edge_flags")
    end
  end
  i, _ = core.edge_flags[_e]
  v = src(_e)
  return sign * weight(core, v, i)
end

function weight(G::AbstractGKMGraph, e::Edge)
  return weight(core(G), e)
end

@doc raw"""
    is2_indep(G::AbstractGKM_graph) -> Bool

Return `true` if `G` is 2-independent, i.e. the weights of every two flags at a vertex are linearly independent.
"""
function is2_indep(G::AbstractGKMGraph)
  return _indep(core(G), 2)
end

@doc raw"""
    is3_indep(G::AbstractGKM_graph) -> Bool

Return `true` if `G` is 3-independent, i.e. the weights of every three flags at a vertex are linearly independent.
# Example
The weights of $\mathbb{P}^3$ at the fixed point $[1:0:0:0]$ are $\{t_i-t_0:i\in\{1, 2, 3\}\}$, which are linearly independent over $\mathbb{C}$.
```jldoctest is3_indep
julia> is3_indep(projective_space(GKM_graph, 3))
true
```
The variety of complete flags in $\mathbb{C}^3$ is an example of a GKM graph that is not 3-independent:
```jldoctest is3_indep
julia> G = flag_variety(GKM_graph, [1, 1, 1])
GKM graph with 6 nodes, valency 3 and axial function:
13 -> 12 => (0, -1, 1)
21 -> 12 => (-1, 1, 0)
23 -> 13 => (-1, 1, 0)
23 -> 21 => (-1, 0, 1)
31 -> 13 => (-1, 0, 1)
31 -> 21 => (0, -1, 1)
32 -> 12 => (-1, 0, 1)
32 -> 23 => (0, -1, 1)
32 -> 31 => (-1, 1, 0)

julia> is3_indep(G)
false
```
!!! warning
    This function throws an error if the valency of `G` is less than 3, since in this case it is not possible to pick three different flags at a vertex.
"""
function is3_indep(G::AbstractGKMGraph)
  return _indep(core(G), 3)
end

function _indep(core::GKMCombinatorialData, k::Int64)

  @req valency(core) >= k "valency is too low"

  val = valency(core)

  for v in 1:n_vertices(core.g)
    # Check all k-tuples of distinct flag indices at vertex v
    for tup in Iterators.product([1:val for _ in 1:k]...)
      # Skip if not strictly increasing (to avoid checking same set multiple times)
      any(i -> tup[i-1] >= tup[i], 2:k) && continue

      # Get the weights of the k flags
      weights = [weight(core, v, tup[i]) for i in 1:k]

      if rank(matrix(weights)) < k
        return false
      end
    end
  end

  return true
end
