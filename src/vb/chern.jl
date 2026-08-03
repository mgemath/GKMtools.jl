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

@doc raw"""
    chern_class(E::AbstractGKMVectorBundle, k::Integer)
    chern_class(G::AbstractGKMGraph, k::Integer)

Return the `k`-th equivariant Chern class of a smooth or orbifold GKM vector
bundle `E`. For a graph `G`, return the Chern class of its tangent bundle.

The result is a localized [`GKMClass`](@ref), whose restriction at each vertex
is the `k`-th elementary symmetric polynomial in the fiber weights. The zeroth
Chern class is the unit, and classes above the bundle rank are zero.

Throw an error when `k` is negative or the torus ranks of the bundle and its
base differ.
"""
function Oscar.chern_class(
  E::AbstractGKMVectorBundle,
  k::Integer,
)
  @req k >= 0 """
  Chern class is only defined for non-negative index
  """

  G = baseof(E)
  basis = gens_cohomRing(G)

  k == 0 && return one(first(basis))
  k > rank(E) && return zero(first(basis))

  parameters = gens_coeffRing(G)
  result = zero(first(basis))

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

@doc raw"""
    first_chern_class(E::AbstractGKMVectorBundle)
    first_chern_class(G::AbstractGKMGraph)

Return the first equivariant Chern class of `E`, or of the tangent bundle of
`G`. This is equivalent to `chern_class(E, 1)` or `chern_class(G, 1)`.
"""
first_chern_class(
  E::AbstractGKMVectorBundle,
) = chern_class(E, 1)

"""
    chern_classes(E::AbstractGKMVectorBundle)
    chern_classes(G::AbstractGKMGraph)

Return all equivariant Chern classes from degree zero through the rank. For a
graph, compute the classes of its tangent bundle. The first element is the
unit class.
"""
function Oscar.chern_classes(
  E::AbstractGKMVectorBundle,
)
  return [
    chern_class(E, k)
    for k in 0:rank(E)
  ]
end

@doc raw"""
    total_chern_class(E::AbstractGKMVectorBundle)
    total_chern_class(G::AbstractGKMGraph)

Return the total equivariant Chern class, defined as the sum of all Chern
classes of `E`. For a graph, compute the total Chern class of its tangent
bundle.
"""
function total_chern_class(
  E::AbstractGKMVectorBundle,
)
  unit = first(gens_cohomRing(baseof(E)))

  return sum(
    chern_classes(E);
    init=zero(unit),
  )
end

function Oscar.chern_class(
  G::AbstractGKMGraph,
  k::Integer,
)
  return chern_class(tangent_bundle(G), k)
end

Oscar.chern_classes(
  G::AbstractGKMGraph,
) = chern_classes(tangent_bundle(G))

total_chern_class(
  G::AbstractGKMGraph,
) = total_chern_class(tangent_bundle(G))

first_chern_class(
  G::AbstractGKMGraph,
) = chern_class(G, 1)
