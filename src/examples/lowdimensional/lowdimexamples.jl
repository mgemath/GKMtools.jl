@doc raw"""
    gkm_2d(w::AbstractMatrix{<:Integer}) -> GKMGraph

Return the 2-valent GKM cyclic connected GKM graph whose vertices are $1,2,\dots,n$ and whose edges
are $(1, 2), (2, 3), ..., (n, 1)$.
The weights of those edges are given by the rows of the matrix `w`.

# Example

The following example is the GKM graph from [GKZ22; Example 2.44, left figure](@cite), which cannot come from a Hamiltonian action.
One way of seeing this is that the combinatorial Betti numbers are not the geometric Betti numbers of any connected space.

```jldoctest
julia> G = gkm_2d([1 0; 0 1; -1 0; 0 -1; 1 0; 0 1; -1 0; 0 -1;])
GKM graph with 8 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (0, -1)
4 -> 3 => (1, 0)
5 -> 4 => (0, 1)
6 -> 5 => (-1, 0)
7 -> 6 => (0, -1)
8 -> 1 => (0, -1)
8 -> 7 => (1, 0)

julia> betti_numbers(G)
3-element Vector{Int64}:
 2
 4
 2
```
"""
function _gkm_graph_from_weighted_edges(n::Int, r::Int, weighted_edges)
  g = Graph{Undirected}(n)
  M = free_module(ZZ, r)
  flags = [FlagWeight{ZZRingElem}[] for _ in 1:n]
  edge_flags = Dict{Edge,Tuple{Int,Int}}()
  for (u, v, coordinates) in weighted_edges
    @req !has_edge(g, u, v) "multiple edges are not supported"
    add_edge!(g, u, v)
    e = Edge(u, v) in edges(g) ? Edge(u, v) : Edge(v, u)
    weight_at_u = M(coordinates)
    weight_at_source = src(e) == u ? weight_at_u : -weight_at_u
    push!(flags[src(e)], FlagWeight{ZZRingElem}(weight_at_source))
    push!(flags[dst(e)], FlagWeight{ZZRingElem}(-weight_at_source))
    edge_flags[e] = (length(flags[src(e)]), length(flags[dst(e)]))
  end
  labels = Vertex.(string.(1:n))
  core = GKMCombinatorialData{ZZRingElem,Vertex,FlagWeight{ZZRingElem}}(g, M, labels, flags, edge_flags)
  return GKMGraph{ZZRingElem,Vertex,FlagWeight{ZZRingElem}}(
    core, build_gkm_connection(core), create_cohomology(r, n), nothing,
  )
end

function _gkm_2d_with_extra_edges(w::AbstractMatrix{<:Integer}, extra_edges)
  n, r = size(w)
  cycle_edges = [(i, i % n + 1, collect(w[i, :])) for i in 1:n]
  return _gkm_graph_from_weighted_edges(n, r, vcat(cycle_edges, extra_edges))
end

function gkm_2d(w::AbstractMatrix{<:Integer})::GKMGraph
  n, r = size(w)
  @req n >= 3 "Need at least three vertices for 2d GKM graph"
  weighted_edges = [(i, i % n + 1, collect(w[i, :])) for i in 1:n]
  return _gkm_graph_from_weighted_edges(n, r, weighted_edges)
end

@doc raw"""
    gkm_3d_positive_non_toric(i::Integer) -> GKMGraph

Return the `i`-th GKM graph from [CK23; Appendix A](@cite)'s classification of 3-valent Hamiltonian positive GKM-graphs
with 2-dimensional torus-action that are not projections of GKM graphs coming from smooth projective polytopes.

The argument `i` runs from 1 ot 7 and is the index in the list.

# Example
We reproduce here the Betti numbers and the integrals $\int_M (c_1(M))^3$ as listed in [CK23; Appendix A](@cite).
```jldoctest
julia> for i in 1:7
           G = gkm_3d_positive_non_toric(i)
           integral = integrate(first_chern_class(G)^3, G)
           println("Graph $i: Betti numbers = $(betti_numbers(G)), [M] . (c_1(M))^3 = $integral")
         end
Graph 1: Betti numbers = [1, 1, 1, 1], [M] . (c_1(M))^3 = 54
Graph 2: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 30
Graph 3: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 40
Graph 4: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 46
Graph 5: Betti numbers = [1, 2, 2, 1], [M] . (c_1(M))^3 = 48
Graph 6: Betti numbers = [1, 3, 3, 1], [M] . (c_1(M))^3 = 26
Graph 7: Betti numbers = [1, 3, 3, 1], [M] . (c_1(M))^3 = 38
```
"""
function gkm_3d_positive_non_toric(i::Integer)::GKMGraph
  if i == 1
    return _gkm_2d_with_extra_edges([1 1; -1 1; -1 -1; 1 -1], [(1, 3, [0, 1]), (4, 2, [1, 0])])
  elseif i == 2
    return _gkm_2d_with_extra_edges([2 0; 0 1; -2 2; -1 0; 0 -2; 1 -1], [(1, 4, [0, 1]), (2, 5, [-1, 1]), (3, 6, [-1, 0])])
  elseif i == 3
    return _gkm_2d_with_extra_edges([2 0; 1 1; -1 1; -2 0; -1 -1; 1 -1], [(1, 5, [0, 1]), (2, 4, [0, 1]), (3, 6, [-1, 0])])
  elseif i == 4
    return _gkm_2d_with_extra_edges([1 1; 0 1; -1 1; -1 0; -1 -1; 1 -1], [(1, 4, [0, 1]), (2, 5, [-1, 1]), (3, 6, [-1, 0])])
  elseif i == 5
    return _gkm_2d_with_extra_edges([1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1], [(1, 4, [0, 1]), (2, 5, [-1, 1]), (3, 6, [-1, 0])])
  elseif i == 6
    return _gkm_2d_with_extra_edges([2 0; 1 1; 0 2; -1 1; -2 0; -1 -1; 0 -2; 1 -1], [(1, 6, [0, 1]), (2, 5, [0, 1]), (3, 8, [-1, 0]), (4, 7, [-1, 0])])
  elseif i == 7
    return _gkm_2d_with_extra_edges([1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1], [(1, 5, [0, 1]), (2, 7, [-1, 1]), (3, 6, [-1, 1]), (4, 8, [-1, 0])])
  else
    @req false "Index must be between 1 and 7."
  end
end

@doc raw"""
    gkm_3d_twisted_flag() -> GKMGraph

Return the 3-valent GKM graph of the twisted flag manifolds of Eschenburg, Tolman, and Woodward
(see [GKZ20; Example 4.8](@cite) and references therein).

# Example
Note that the resulting GKM graph does not occur in the output of `gkm_3d_positive_non_toric()` since one edge has
non-positive Chern number.

```jldoctest
julia> G = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)

julia> print_curve_classes(G)
2 -> 1: (0, 1), Chern number: 4
3 -> 2: (-1, 1), Chern number: 2
4 -> 1: (1, 0), Chern number: 2
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 2
5 -> 4: (-1, 1), Chern number: 2
6 -> 1: (1, 1), Chern number: 6
6 -> 3: (1, 0), Chern number: 2
6 -> 5: (0, 1), Chern number: 4
```
"""
function gkm_3d_twisted_flag()::GKMGraph
  return _gkm_2d_with_extra_edges([0 1; -1 0; 1 -1; 0 1; -1 0; 1 -1], [(1, 4, [-1, 2]), (2, 5, [-1, 1]), (3, 6, [-2, 1])])
end

#TODO: document this, referring to [GKZ20, Prop. 4.5]'s classification of signed 3D GKM fibrations. 
function gkm_3d_fibration(w::AbstractMatrix{<:Integer}, k::AbstractVector{<:Integer}, twisted::Bool)::GKMGraph
  n, r = size(w)
  @req n >= 3 "Base 2D GKM graph needs at least 3 vertices"
  @req length(k) == n "there must be one integer k for each base vertex"
  weighted_edges = Tuple[]
  for i in 1:(n - 1)
    push!(weighted_edges, (i, i + 1, collect(w[i, :])))
    push!(weighted_edges, (i + n, i + n + 1, collect(w[i, :])))
  end
  for i in 2:n
    gi = k[i] .* w[i - 1, :] .- k[i - 1] .* w[i, :]
    push!(weighted_edges, (i, i + n, collect(gi)))
  end
  if twisted
    push!(weighted_edges, (n, n + 1, collect(w[n, :])))
    push!(weighted_edges, (2 * n, 1, collect(w[n, :])))
    gi = k[1] .* w[n, :] .+ k[n] .* w[1, :]
  else
    push!(weighted_edges, (n, 1, collect(w[n, :])))
    push!(weighted_edges, (2 * n, n + 1, collect(w[n, :])))
    gi = k[1] .* w[n, :] .- k[n] .* w[1, :]
  end
  push!(weighted_edges, (1, n + 1, collect(gi)))
  return _gkm_graph_from_weighted_edges(2 * n, r, weighted_edges)
end
