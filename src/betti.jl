@doc raw"""
    betti_numbers(G::AbstractGKM_graph) -> Vector{Int64}

Return the array `betti_numbers` such that `betti_numbers[i+1]` is the 2i-th combinatorial Betti number for i from 0 to the valency of `G`, as defined in [GZ01; section 1.3](@cite).
!!! note
    * `i` ranges from 0 to the valency of `G`, that can be obtained by `valency(G)`.
    * From [GZ01; Theorem 1.3.2](@cite), the combinatorial Betti numbers equal the Betti numbers of the underlying GKM space if the torus action is Hamiltonian.
      This holds automatically for smooth projective varieties with algebraic torus action (cf. [GFK12; Example 8.1 (ii)](@cite)).

!!! warning
    * `betti_numbers[1]` is the 0-th Betti number, since Julia arrays are 1-based and not 0-based.
    * Currently only implemented for compact GKM spaces.

# Examples

```jldoctest
julia> H6 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 6));

julia> betti_numbers(H6)
3-element Vector{Int64}:
 1
 2
 1
```

""" 
function Oscar.betti_numbers(G::AbstractGKM_graph)::Vector{Int64}

  @req is_compact(G) "betti_numbers only defined for compact GKM spaces"

  # weight type parent ring
  PR = parent(zero(G.weightType))
  # weight type
  WT = G.weightType

  for counter in 1:10^8 # arbitrary maximum number of attempts to avoid "while true"

    xi = PR.(rand(Int, rank_torus(G)))  #TODO: find something without using random numbers
    wxi = Dict{Edge, WT}() # wxi stands for weight[e](xi)
    isPolarizing = true

    # calculate weight[e](xi) for all edges e
    for e in edges(G.g)

      wxi[e] = wxi[reverse(e)] = 0

      for j in 1:rank_torus(G)
        wxi[e] += xi[j] * _w(G, e)[j]
        wxi[reverse(e)] += xi[j] * _w(G, reverse(e))[j]
      end

      if wxi[e] == 0 || wxi[reverse(e)] == 0
        isPolarizing = false
        break
      end
    end
    (!isPolarizing) && continue

    # from here on, xi is known to be polarizing.
    #betti = Dict{Int, Int}([i => 0 for i in 0:valency(G)])
    betti = zeros(Int, valency(G)+1)

    for v in 1:n_vertices(G.g)

      i = count(w -> wxi[Edge(v,w)] < 0, all_neighbors(G.g, v) )
      betti[i+1] += 1
    end
    return betti
  end

end



@doc raw"""
    index_periodic_betti(G::AbstractGKM_graph) -> Vector{Int64}

Return the index-periodic Betti numbers of a compact GKM graph.

For a GKM graph with Fano index $p$, this function computes the sums of Betti numbers
grouped by their residue class modulo $p$. Specifically, the $i$-th entry (for $i$ from $1$ to $p$)
contains the sum of all Betti numbers $b_j$ where $j \equiv i$ (mod $p$).

This is as defined in [belmans2025adediagramshodgetatehyperplane; before Theorem 2.2](@cite).

!!! note
    * The returned vector has length equal to the Fano index $p$.
    * Entry i corresponds to the sum of Betti numbers with index congruent $i$ (mod $p$).
    * Requires the GKM graph to be compact.

# Examples
```jldoctest index_periodic_betti
julia> P2 = projective_space(GKM_graph, 2);

julia> index_periodic_betti(P2)
3-element Vector{Int64}:
 1
 1
 1
```
"""
function index_periodic_betti(G::AbstractGKM_graph)::Vector{Int64}
  @req is_compact(G) "index_periodic_betti only defined for compact GKM spaces"

  p = Int64(fano_index(G))
  @req !iszero(p) "Index periodic betti numbers not defined for Calabi Yau spaces (Fano index zero)."
  betti = betti_numbers(G)

  # Initialize the index-periodic Betti numbers
  periodic_betti = zeros(Int64, p)

  # Sum Betti numbers by residue class modulo p
  for (idx, b) in enumerate(betti)
    # idx-1 gives the actual Betti number index (0-based)
    # We want to group by (idx-1) mod p
    residue_class = mod(idx - 2, p) + 1  # +1 for indexing from 1 to p
    periodic_betti[residue_class] += b
  end

  return periodic_betti
end

@doc raw"""
    QH_ss_check_GLLXBR(G::AbstractGKM_graph)::Bool

Return true if $G$ satisfies conditions (1) and (2) in [belmans2025adediagramshodgetatehyperplane; Theorem 2.2](@cite).
"""
function QH_ss_check_GLLXBR(G::AbstractGKM_graph)::Bool
  b = index_periodic_betti(G)
  p = length(b)
  for i in 1:p
    for d in 1:p
      if b[mod(i*d - 1, p) + 1] < b[i]
        return false
      end
    end
  end
  return true
end