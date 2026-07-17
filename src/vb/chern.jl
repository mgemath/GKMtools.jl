function _bundle_weight_class(
  E::AbstractGKMVectorBundle,
  v::Int,
  i::Int,
  t::Vector{T},
) where {T}
  G = baseof(E)

  @req rank(E.M) == rank_torus(G) """
  Chern classes are currently supported when bundle and base have the same torus rank
  """

  res = zero(t[1])
  w = fiber_weight(E, v, i)

  for j in 1:rank_torus(G)
    res += w[j] * t[j]
  end

  return res
end

function _elementary_symmetric_class(values, k::Int)
  @req k >= 0 """
  Chern class is only defined for non-negative index
  """

  k == 0 && return one(first(values))
  k > length(values) && return zero(first(values))

  res = zero(first(values))

  for combination in Combinatorics.combinations(
    eachindex(values),
    k,
  )
    res += prod(
      (values[i] for i in combination);
      init=one(first(values)),
    )
  end

  return res
end

"""
    chern_class(E::AbstractGKMVectorBundle, k::Integer)

Return the `k`-th equivariant Chern class of a smooth or orbifold GKM
vector bundle in the localized cohomology ring of its base.
"""
function Oscar.chern_class(
  E::AbstractGKMVectorBundle,
  k::Integer,
)
  @req k >= 0 """
  Chern class is only defined for non-negative index
  """

  G = baseof(E)
  H = get_cohomology(G)

  k == 0 && return one(H.localized_cohomology)
  k > rank(E) && return zero(H.localized_cohomology)

  basis = gens_cohomRing(H)
  parameters = gens_coeffRing(H.localized_cohomology)
  result = zero(H.localized_cohomology)

  for v in 1:num_vertices(G)
    local_weights = [
      _bundle_weight_class(E, v, i, parameters)
      for i in 1:rank(E)
    ]

    result +=
      _elementary_symmetric_class(
        local_weights,
        Int(k),
      ) * basis[v]
  end

  return result
end

first_chern_class(
  E::AbstractGKMVectorBundle,
) = chern_class(E, 1)

function Oscar.chern_classes(
  E::AbstractGKMVectorBundle,
)
  return [
    chern_class(E, k)
    for k in 0:rank(E)
  ]
end

function total_chern_class(
  E::AbstractGKMVectorBundle,
)
  H = get_cohomology(baseof(E))

  return sum(
    chern_classes(E);
    init=zero(H.localized_cohomology),
  )
end