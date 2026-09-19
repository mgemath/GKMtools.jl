@doc raw"""
    vector_bundle_O(X::WeightedProjectiveSpace, degrees; small_torus=false)

Return the direct sum `O(degrees[1]) ⊕ ⋯ ⊕ O(degrees[end])` on `X`.
"""
function vector_bundle_O(
  X::WeightedProjectiveSpace,
  degrees::AbstractVector{<:Integer};
  small_torus::Bool=false,
)
  @req !isempty(degrees) "Need at least one line-bundle degree"
  return direct_sum(
    (line_bundle_O(X, degree; small_torus=small_torus) for degree in degrees)...,
  )
end

@doc raw"""
    vector_bundle_O(n::Integer, degrees; small_torus=false, enlarge_torus=false)

Return a direct sum of line bundles on ordinary projective `n`-space.

If `enlarge_torus` is `true`, add one independent fibre-scaling character for
each summand. This recovers the linearization used by the former implementation.
"""
function vector_bundle_O(
  n::Integer,
  degrees::AbstractVector{<:Integer};
  small_torus::Bool=false,
  enlarge_torus::Bool=false,
)
  @req n >= 1 "The dimension must be positive"
  @req !isempty(degrees) "Need at least one line-bundle degree"
  @req !small_torus "vector_bundle_O(n, degrees) currently uses the coordinate torus; use small_torus=false"

  G = projective_space(GKMGraph, Int(n))
  if enlarge_torus
    old_lattice = lattice(G)
    M = free_module(base_ring(old_lattice), rank(old_lattice) + length(degrees))
    inclusion = hom(old_lattice, M, gens(M)[1:rank(old_lattice)])
    old_core = core(G)
    new_flags = [
      [FlagWeight(inclusion(weight(G, v, i))) for i in 1:valency(G)]
      for v in vertices(G)
    ]
    G = gkm_graph(GKMCombinatorialData(
      graph(G), M, old_core.labels, new_flags, old_core.edge_flags,
    ))

    basis = gens(M)
    weights = Matrix{eltype(basis)}(undef, num_vertices(G), length(degrees))
    for (j, degree) in enumerate(degrees)
      weights[1, j] = basis[rank(old_lattice) + j]
      for v in 2:num_vertices(G)
        weights[v, j] = weights[1, j] - degree * weight(G, Edge(1, v))
      end
    end
    return vector_bundle(G, M, hom(M, M, basis), weights)
  end

  return direct_sum(
    (_line_bundle_O(G, degree) for degree in degrees)...,
  )
end

@doc raw"""
    vector_bundle_O(n::AbstractVector{<:Integer}, degrees; small_torus=false, enlarge_torus=false)

Return a direct sum of line bundles on `P^n[1] × ⋯ × P^n[end]`.
Each entry of `degrees` is an integer vector of length `length(n)`, giving
one summand's degree on each projective-space factor. For example,
`vector_bundle_O([1, 2], [[2, 3], [-1, 0]])` constructs
`O(2, 3) ⊕ O(-1, 0)` on `P¹ × P²`.

The coordinate torus is used. With `enlarge_torus=true`, add one independent
fibre-scaling character per summand. Otherwise the linearization is the
sum of the factor linearizations used by `vector_bundle_O(::Integer, ...)`.
"""
function vector_bundle_O(
  n::AbstractVector{<:Integer},
  degrees::AbstractVector{<:AbstractVector{<:Integer}};
  small_torus::Bool=false,
  enlarge_torus::Bool=false,
)
  @req !isempty(n) "Need at least one projective-space factor"
  @req all(d -> d >= 1, n) "The dimensions must be positive"
  @req !isempty(degrees) "Need at least one line-bundle multidegree"
  @req all(d -> length(d) == length(n), degrees) "Each multidegree must have one entry per projective-space factor"
  @req !small_torus "vector_bundle_O(n, degrees) currently uses the coordinate torus; use small_torus=false"

  factors = [projective_space(GKMGraph, Int(d)) for d in n]
  G = foldl(*, factors)
  old_lattice = lattice(G)
  M = enlarge_torus ?
    free_module(base_ring(old_lattice), rank(old_lattice) + length(degrees)) :
    old_lattice
  basis = gens(M)
  if enlarge_torus
    inclusion = hom(old_lattice, M, basis[1:rank(old_lattice)])
    old_core = core(G)
    new_flags = [
      [FlagWeight(inclusion(weight(G, v, i))) for i in 1:valency(G)]
      for v in vertices(G)
    ]
    G = gkm_graph(GKMCombinatorialData(
      graph(G), M, old_core.labels, new_flags, old_core.edge_flags,
    ))
  end

  # Product vertices vary fastest in the first factor; lattice blocks follow
  # the same factor order.
  offsets = cumsum(vcat(0, [rank(lattice(F)) for F in factors]))
  inclusions = [
    hom(lattice(F), M, basis[offsets[i] + 1:offsets[i + 1]])
    for (i, F) in enumerate(factors)
  ]
  origin = sum(basis[offsets[i] + 1] for i in eachindex(factors))
  weights = Matrix{eltype(basis)}(undef, num_vertices(G), length(degrees))
  for (v, point) in enumerate(Iterators.product((vertices(F) for F in factors)...))
    for (j, degree) in enumerate(degrees)
      w = enlarge_torus ? basis[rank(old_lattice) + j] : origin
      for (i, F) in enumerate(factors)
        point[i] == 1 && continue
        w -= degree[i] * inclusions[i](weight(F, Edge(1, point[i])))
      end
      weights[v, j] = w
    end
  end
  return vector_bundle(G, M, hom(M, M, basis), weights)
end

gkm_vector_bundle_of_toric(L::ToricLineBundle) = gkm_vector_bundle_of_toric([L])

@doc raw"""
    gkm_vector_bundle_of_toric(E::AbstractVector{<:ToricLineBundle}) -> GKMVectorBundle

Return the GKM vector bundle associated with the direct sum of the toric line
bundles in `E`. All bundles must be defined on the same smooth projective
toric variety.

# Example
```jldoctest
julia> X = hirzebruch_surface(NormalToricVariety, 2);

julia> E = [toric_line_bundle(X, [1, 0]), toric_line_bundle(X, [0, 1])];

julia> V = gkm_vector_bundle_of_toric(E);

julia> (rank(V), num_vertices(baseof(V)))
(2, 4)
```
"""
function gkm_vector_bundle_of_toric(E::AbstractVector{<:ToricLineBundle})
  @req !isempty(E) "Need at least one toric line bundle"
  X = toric_variety(first(E))
  @req all(L -> toric_variety(L) == X, E) "toric line bundles must have the same base"
  @req is_projective(X) "toric variety must be projective"
  @req is_smooth(X) "toric variety must be smooth"

  total = total_space(E...)
  total_graph = gkm_graph_of_toric(total)
  G = subgraph(subgraph_from_vertices(total_graph, collect(vertices(total_graph))))
  M = lattice(G)
  weights = Matrix{eltype(gens(M))}(undef, num_vertices(G), length(E))
  basis = gens(M)
  fibre_rank = length(E)

  for j in 1:fibre_rank
    fibre_ray = ray_vector([i == dim(total) - fibre_rank + j for i in 1:dim(total)])
    for sigma_index in 1:num_vertices(G)
      dual_ray = nothing
      for candidate in rays(polarize(maximal_cones(total)[sigma_index]))
        iszero(dot(fibre_ray, candidate)) && continue
        dual_ray = lcm(denominator.(candidate)) * candidate
        break
      end
      @assert !isnothing(dual_ray) "could not determine a fibre weight"
      weights[sigma_index, j] = -sum(
        k -> Int(dot(rays(total)[k], dual_ray)) * basis[k],
        1:n_rays(total);
        init=zero(M),
      )
    end
  end

  GMtoM = hom(M, M, basis)
  return vector_bundle(G, M, GMtoM, weights)
end
