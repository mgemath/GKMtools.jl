@doc raw"""
    betti_numbers(G::AbstractGKMGraph) -> Vector{Int64}

Return the array `betti_numbers` such that `betti_numbers[i+1]` is the
``2i``-th combinatorial Betti number for ``i`` from ``0`` to the valency
of `G`, as defined in [GZ01; section 1.3](@cite).

!!! note
    * `i` ranges from `0` to `valency(G)`.
    * From [GZ01; Theorem 1.3.2](@cite), the combinatorial Betti numbers
      equal the Betti numbers of the underlying GKM space if the torus
      action is Hamiltonian.

!!! warning
    * `betti_numbers[1]` is the 0-th Betti number, since Julia arrays are
      1-based.
    * Currently only implemented for compact GKM spaces.
"""
function Oscar.betti_numbers(G::AbstractGKMGraph)::Vector{Int64}
  @req is_compact(G) "betti_numbers only defined for compact GKM spaces"

  for _ in 1:10^8
    xi = rand(-100:100, rank_torus(G))
    all(iszero, xi) && continue

    wxi = Dict{Tuple{Int,Int},Any}()
    is_polarizing = true

    for e in edges(G)
      forward = _edge_weight_pairing(G, e, xi)
      backward = _edge_weight_pairing(G, Edge(dst(e), src(e)), xi)
      wxi[(src(e), dst(e))] = forward
      wxi[(dst(e), src(e))] = backward

      if iszero(forward) || iszero(backward)
        is_polarizing = false
        break
      end
    end

    is_polarizing || continue

    betti = zeros(Int, valency(G) + 1)
    for v in vertices(G)
      index = count(w -> wxi[(v, w)] < 0, all_neighbors(graph(G), v))
      betti[index + 1] += 1
    end
    return betti
  end

  error("Could not find a polarizing vector for the GKM graph")
end

function _edge_weight_pairing(G::AbstractGKMGraph, e::Edge, xi)
  w = weight(G, e)
  res = zero(w[1])
  for j in 1:rank_torus(G)
    res += xi[j] * w[j]
  end
  return res
end

function _first_chern_weight(G::AbstractGKMGraph, v::Int)
  @req !isempty(flags(G, v)) "Cannot compute Chern weight at a vertex with no flags"
  res = zero(lattice(G))
  for i in 1:length(flags(G, v))
    res += weight(G, v, i)
  end
  return res
end

function _proportional_scalar(numerator_weight, denominator_weight)
  scalar = nothing
  for j in 1:rank(parent(denominator_weight))
    if iszero(denominator_weight[j])
      @req iszero(numerator_weight[j]) "Weights are not proportional"
    else
      q = numerator_weight[j] / denominator_weight[j]
      if isnothing(scalar)
        scalar = q
      else
        @req scalar == q "Weights are not proportional"
      end
    end
  end
  @req !isnothing(scalar) "Cannot divide by the zero edge weight"
  return scalar
end

function _integer_scalar(q)::ZZRingElem
  @req denominator(q) == 1 "1st Chern number is not integral"
  return ZZ(numerator(q))
end

@doc raw"""
    chern_number(e::Edge, G::AbstractGKMGraph) -> ZZRingElem

Return the first Chern number of the invariant curve represented by `e`.
"""
function Oscar.chern_number(e::Edge, G::AbstractGKMGraph)::ZZRingElem
  src_weight = _first_chern_weight(G, src(e))
  dst_weight = _first_chern_weight(G, dst(e))
  scalar = _proportional_scalar(src_weight - dst_weight, weight(G, e))
  return _integer_scalar(scalar)
end

@doc raw"""
    fano_index(G::AbstractGKMGraph) -> ZZRingElem

Return the Fano index of the GKM graph, which is the greatest common divisor of the first Chern numbers of all its edges.

# Examples
```jldoctest fano_index
julia> P2 = projective_space(GKMGraph, 2);

julia> fano_index(P2)
3

julia> F3 = flag_variety(GKMGraph, [1, 1, 1]);

julia> fano_index(F3)
2
```
"""
function fano_index(G::AbstractGKMGraph)::ZZRingElem
  chern_numbers = [abs(chern_number(e, G)) for e in edges(G)]
  isempty(chern_numbers) && return ZZ(0)
  return foldl(gcd, chern_numbers)
end

@doc raw"""
    index_periodic_betti(G::AbstractGKMGraph) -> Vector{Int64}

Return the index-periodic Betti numbers of a compact GKM graph.

For a GKM graph with Fano index ``p``, this function computes the sums of
Betti numbers grouped by their residue class modulo ``p``.
"""
function index_periodic_betti(G::AbstractGKMGraph)::Vector{Int64}
  @req is_compact(G) "index_periodic_betti only defined for compact GKM spaces"

  p = Int64(fano_index(G))
  @req !iszero(p) "Index periodic betti numbers not defined for Calabi Yau spaces (Fano index zero)."
  betti = betti_numbers(G)

  periodic_betti = zeros(Int64, p)
  for (idx, b) in enumerate(betti)
    residue_class = mod(idx - 2, p) + 1
    periodic_betti[residue_class] += b
  end

  return periodic_betti
end


@doc raw"""
    pseudo_index(G::AbstractGKMGraph) -> ZZRingElem

Return the pseudo index of the GKM graph, which is the minimum of the first Chern numbers of all its edges.

# Examples
The examples below show that the pseudo-index differs from the Fano index in general.
```jldoctest pseudo_index
julia> P2 = projective_space(GKMGraph, 2);

julia> fano_index(P2), pseudo_index(P2)
(3, 3)

julia> P3 = projective_space(GKMGraph, 3);

julia> fano_index(P3), pseudo_index(P3)
(4, 4)

julia> G = P2 * P3;

julia> fano_index(G), pseudo_index(G)
(1, 3)

julia> T = gkm_3d_twisted_flag();

julia> fano_index(T), pseudo_index(T)
(2, 0)
```
"""
function pseudo_index(G::AbstractGKMGraph)::ZZRingElem
  chern_nums = [chern_number(e, G) for e in edges(G.g)]
  return minimum(chern_nums)
end


@doc raw"""
    QH_ss_check_GLLXBR(G::AbstractGKMGraph)::Bool

Return true if `G` satisfies conditions (1) and (2) in
[belmans2025adediagramshodgetatehyperplane; Theorem 2.2](@cite).
"""
function QH_ss_check_GLLXBR(G::AbstractGKMGraph)::Bool
  b = index_periodic_betti(G)
  p = length(b)
  for i in 1:p
    for d in 1:p
      if b[mod(i * d - 1, p) + 1] < b[i]
        return false
      end
    end
  end
  return true
end

@doc raw"""
    is_strictly_nef(G::AbstractGKMGraph) -> Bool

Return `true` if and only if the Chern numbers of all curve classes corresponding to
edges of the GKM graph are strictly positive.

# Examples
```jldoctest is_strictly_nef
julia> F3 = flag_variety(GKMGraph, [1,1,1]);

julia> print_curve_classes(F3)
13 -> 12: (0, 1), Chern number: 2
21 -> 12: (1, 0), Chern number: 2
23 -> 13: (1, 1), Chern number: 4
23 -> 21: (0, 1), Chern number: 2
31 -> 13: (1, 0), Chern number: 2
31 -> 21: (1, 1), Chern number: 4
32 -> 12: (1, 1), Chern number: 4
32 -> 23: (1, 0), Chern number: 2
32 -> 31: (0, 1), Chern number: 2

julia> is_strictly_nef(F3)
true

julia> H5 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5))
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (5, 1, 0, -1)
4 -> 1 => (0, 1, 5, -1)
4 -> 3 => (-1, 0, 1, 0)

julia> print_curve_classes(H5)
2 -> 1: (-5, 1), Chern number: -3
3 -> 2: (1, 0), Chern number: 2
4 -> 1: (1, 0), Chern number: 2
4 -> 3: (0, 1), Chern number: 7

julia> is_strictly_nef(H5)
false
```
"""
function is_strictly_nef(G::AbstracGraph)::Bool
  for e in edges(G.g)
    if chern_number(e, G) <= 0
      return false
    end
  end
  return true
end
