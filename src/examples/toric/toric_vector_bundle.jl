"""
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

"""
    vector_bundle_O(n::Integer, degrees; small_torus=false)

Return a direct sum of line bundles on ordinary projective `n`-space.
"""
function vector_bundle_O(
  n::Integer,
  degrees::AbstractVector{<:Integer};
  small_torus::Bool=false,
)
  @req n >= 1 "The dimension must be positive"
  @req !isempty(degrees) "Need at least one line-bundle degree"
  @req !small_torus "vector_bundle_O(n, degrees) currently uses the coordinate torus; use small_torus=false"

  G = projective_space(GKMGraph, Int(n))
  return direct_sum(
    (_line_bundle_O(G, degree) for degree in degrees)...,
  )
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
