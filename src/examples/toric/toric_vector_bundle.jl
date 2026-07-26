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
