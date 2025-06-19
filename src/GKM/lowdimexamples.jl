@doc raw"""
    gkm_2d(w::Matrix{Int64}) -> AbstractGKM_graph

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
function gkm_2d(w::Matrix{Int64})::AbstractGKM_graph
  s = size(w)
  n = s[1]
  @req n >= 3 "Need at least three vertices for 2d GKM graph"
  r = s[2]
  labels = ["$i" for i in 1:n]
  G = empty_gkm_graph(n, r, labels)

  for i in 1:n
    add_edge!(G, i, i % n + 1, G.M(w[i,:]))
  end
  return G
end

@doc raw"""
    gkm_3d_positive_non_toric(i::Int64) -> AbstractGKM_graph

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
function gkm_3d_positive_non_toric(i::Int64)::AbstractGKM_graph
  if i == 1
    G = gkm_2d([1 1; -1 1; -1 -1; 1 -1])
    add_edge!(G, 1, 3, G.M([0, 1]))
    add_edge!(G, 4, 2, G.M([1, 0]))
    return G
  elseif i == 2
    G = gkm_2d([2 0; 0 1; -2 2; -1 0; 0 -2; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 3
    G = gkm_2d([2 0; 1 1; -1 1; -2 0; -1 -1; 1 -1])
    add_edge!(G, 1, 5, G.M([0, 1]))
    add_edge!(G, 2, 4, G.M([0, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 4
    G = gkm_2d([1 1; 0 1; -1 1; -1 0; -1 -1; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 5
    G = gkm_2d([1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1])
    add_edge!(G, 1, 4, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 0]))
    return G
  elseif i == 6
    G = gkm_2d([2 0; 1 1; 0 2; -1 1; -2 0; -1 -1; 0 -2; 1 -1])
    add_edge!(G, 1, 6, G.M([0, 1]))
    add_edge!(G, 2, 5, G.M([0, 1]))
    add_edge!(G, 3, 8, G.M([-1, 0]))
    add_edge!(G, 4, 7, G.M([-1, 0]))
    return G
  elseif i == 7
    G = gkm_2d([1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1])
    add_edge!(G, 1, 5, G.M([0, 1]))
    add_edge!(G, 2, 7, G.M([-1, 1]))
    add_edge!(G, 3, 6, G.M([-1, 1]))
    add_edge!(G, 4, 8, G.M([-1, 0]))
    return G
  else
    @req false "Index must be between 1 and 7."
  end
end

@doc raw"""
    gkm_3d_twisted_flag() -> AbstractGKM_graph

Return the 3-valent GKM graph of the twisted flag varieties of Eschenburg, Tolman, and Woodward
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
function gkm_3d_twisted_flag()::AbstractGKM_graph
  G = gkm_2d([0 1; -1 0; 1 -1; 0 1; -1 0; 1 -1])
  add_edge!(G, 1, 4, G.M([-1, 2]))
  add_edge!(G, 2, 5, G.M([-1, 1]))
  add_edge!(G, 3, 6, G.M([-2, 1]))
  return G
end

#TODO: document this, referring to [GKZ20, Prop. 4.5]'s classification of signed 3D GKM fibrations. 
function gkm_3d_fibration(w::Matrix{Int64}, k::Vector{Int64}, twisted::Bool)::AbstractGKM_graph
  s = size(w)
  n = s[1]
  @req n >= 3 "Base 2D GKM graph needs at least 3 vertices"
  r = s[2]
  labels = ["$i" for i in 1:2*n]
  G = empty_gkm_graph(2*n, r, labels)

  # add most horizontal edges
  for i in 1:(n-1)
    a = i
    b = i % n + 1
    add_edge!(G, a, b, G.M(w[i,:]))
    add_edge!(G, a + n, b + n, G.M(w[i,:]))
  end
  
  # add most vertical edges
  for i in 2:n
    gi = k[i]*G.M(w[i-1,:]) - k[i-1]*G.M(w[i,:])
    add_edge!(G, i, i+n, gi)
  end

  if twisted
    # Horizontal edges are 1-2-3-...-n-(n+1)-...-(2n)-1
    add_edge!(G, n, n + 1, G.M(w[n,:]))
    add_edge!(G, 2*n, 1, G.M(w[n,:]))

    gi = k[1]*G.M(w[n,:]) + k[n]*G.M(w[1,:])
    add_edge!(G, 1, n+1, gi)
  else
    # Horizontal edges are 1-2-3-...-n-1 & (n+1)-(n+2)-...-(2n)-(n+1)
    add_edge!(G, n, 1, G.M(w[n,:]))
    add_edge!(G, 2*n, n+1, G.M(w[n,:]))

    gi = k[1]*G.M(w[n,:]) - k[n]*G.M(w[1,:])
    add_edge!(G, 1, n+1, gi)
  end

  return G
end