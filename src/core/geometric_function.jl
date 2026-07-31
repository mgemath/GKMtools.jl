
rank_torus(core::GKMCombinatorialData) = rank(core.M)

@doc raw"""
    rank_torus(G::AbstractGKMGraph) -> Int64

Return the rank of the torus acting on `G`. That is, the rank of the character group.

# Examples
By default, the torus acting on $\mathbb{P}^n$ is $(\mathbb{C}^\times)^{n+1}$, acting by rescaling the homogeneous coordinates.
```jldoctest rank_torus
julia> P3 = projective_space(GKMGraph, 3);

julia> rank_torus(P3)
4
```
Taking products adds the rank:
```jldoctest rank_torus
julia> H6 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 6));

julia> rank_torus(H6)
4
julia> rank_torus(H6 * P3)
8
```
"""
rank_torus(G::AbstractGKMGraph) = rank_torus(core(G))

@doc raw"""
    valency(G::AbstractGKMGraph) -> Int64

Return the valency of `G`, i.e. the number of flags at each vertex.

!!! warning
    This function does not check if `G` is a valid GKM graph (use `isvalid` to check this).
    In particular, it does not check if every vertex has the same degree.
    The returned value is the degree of vertex `1`.

# Example:
The valency of the GKM graph of $\mathbb{P}^3$ is 3, since all of the fixed points $[1:0:0:0], \dots, [0:0:0:1]$ are connected to each other
via some $T$-invariant $\mathbb{P}^1$'s. For example, $[1:0:0:0]$ and $[0:1:0:0]$ are connected by $\{[x:y:0:0] : x,y\in\mathbb{C}\}$.
```jldoctest valency
julia> valency(projective_space(GKMGraph, 3))
3
julia> valency(grassmannian(GKMGraph, 2, 4)) # The Grassmannian of 2-planes in C^4
4
julia> valency(flag_variety(GKMGraph, [1, 1, 1, 1])) # The variety of full flags in C^4
6
```
"""
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

@doc raw"""
    is_compact(G::AbstractGKMGraph) -> Bool

Return `true` if `G` is compact, i.e. all flags at all vertices are associated with edges (no standalone flags).

# Example
```jldoctest is_compact
julia> G = projective_space(GKMGraph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> is_compact(G)
true

julia> A3 = gkm_graph_of_toric(affine_space(NormalToricVariety, 3)) # affine space
GKM graph with 1 nodes, valency 3 and axial function:
Standalone flags:
Algorithmic connection for GKM graph with 1 nodes and valency 3

julia> is_compact(A3)
false
```
"""
function is_compact(G::AbstractGKMGraph)
  val = valency(G)
  for v in vertices(G)
    if length(neighbors(graph(G), v)) != val
      return false
    end
  end
  return true
end

@doc raw"""
    is2_indep(G::AbstractGKMGraph) -> Bool

Return `true` if `G` is 2-independent, i.e. the weights of every two flags at a vertex are linearly independent.
"""
function is2_indep(G::AbstractGKMGraph)
  return _indep(core(G), 2)
end

@doc raw"""
    is3_indep(G::AbstractGKMGraph) -> Bool

Return `true` if `G` is 3-independent, i.e. the weights of every three flags at a vertex are linearly independent.
# Example
The weights of $\mathbb{P}^3$ at the fixed point $[1:0:0:0]$ are $\{t_i-t_0:i\in\{1, 2, 3\}\}$, which are linearly independent over $\mathbb{C}$.
```jldoctest is3_indep
julia> is3_indep(projective_space(GKMGraph, 3))
true
```
The variety of complete flags in $\mathbb{C}^3$ is an example of a GKM graph that is not 3-independent:
```jldoctest is3_indep
julia> G = flag_variety(GKMGraph, [1, 1, 1])
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

@doc raw"""
    gkm_independence(G::AbstractGKMGraph) -> Int64

Return the maximum integer in `i` in `0,1,...,valency(G)` such that the GKM graph $G$ is $i$-independent.
The GKM graph is $i$-independent if for each vertex $v$, each $i$-tuple of flags at $v$ has linearly independent axial function values.

# Example
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

julia> gkm_independence(G)
2
```
"""
function gkm_independence(G::GKMCombinatorialData)
  val = valency(G)
  max_indep = 0
  for i in 1:val
    if _indep(G, i)
      max_indep = i
    else
      break
    end
  end
  return max_indep
end

gkm_independence(G::AbstractGKMGraph) = gkm_independence(core(G))