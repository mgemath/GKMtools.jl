@doc raw"""
    schubert_basis(G::AbstractGKMGraph; representation=:localized) -> Vector

Return Billey's equivariant Schubert basis of a homogeneous GKM variety
`G = G/P`. The result is grouped by complex codimension: `result[d + 1]`
contains the classes of codimension `d` (cohomological degree `2d`). Within a
group, classes follow the fixed-point ordering of `G`.

Set `representation=:polynomial` to return [`GKMSpline`](@ref) objects over
`H_T^*(pt)`. The default `:localized` preserves the historical fraction-field
representation.
"""
function schubert_basis(
  G::AbstractGKMGraph{R,V,F};
  representation::Symbol=:localized,
) where {R,V<:GeneralizedFlagVertex,F}
  _check_schubert_representation(representation)
  dimensions = length.(flag.(labels(G)))
  classes = [schubert_basis(G, v; representation) for v in 1:num_vertices(G)]
  basis = [eltype(classes)[] for _ in 0:maximum(dimensions)]
  for vertex in eachindex(classes)
    push!(basis[dimensions[vertex] + 1], classes[vertex])
  end
  return basis
end

function schubert_basis(::AbstractGKMGraph; representation::Symbol=:localized)
  _check_schubert_representation(representation)
  throw(ArgumentError("Billey's formula requires a homogeneous GKM graph with Weyl-group vertices"))
end

function _check_schubert_representation(representation::Symbol)
  representation in (:localized, :polynomial) || throw(ArgumentError(
    "representation must be :localized or :polynomial",
  ))
  return representation
end

"""Alias for [`schubert_basis`](@ref)."""
billey_schubert_basis(G::AbstractGKMGraph; kwargs...) = schubert_basis(G; kwargs...)

function schubert_basis(
  G::AbstractGKMGraph{R,V,F},
  vertex_label::String;
  representation::Symbol=:localized,
) where {R,V<:GeneralizedFlagVertex,F}
  return schubert_basis(G, find_vertex_index(vertex_label, G); representation)
end

function schubert_basis(
  G::AbstractGKMGraph{R,V,F},
  vertex::Integer;
  representation::Symbol=:localized,
) where {R,V<:GeneralizedFlagVertex,F}
  _check_schubert_representation(representation)
  1 <= vertex <= num_vertices(G) || throw(BoundsError(labels(G), vertex))

  u = flag(labels(G)[vertex])
  t = collect(gens(equivariant_coefficient_ring(G)))
  values = [
    _billey_localization(G, u, flag(labels(G)[v]), t)
    for v in 1:num_vertices(G)
  ]
  polynomial = polynomial_class(G, values)
  return representation == :polynomial ? polynomial : localize(polynomial)
end

# Dynamic programming collects the root products of reduced subwords, avoiding
# an explicit enumeration of all 2^length(v) subwords.
function _billey_localization(G::AbstractGKMGraph, u, v, t=gens_coeffRing(G))
  W = parent(v)
  W == parent(u) || throw(ArgumentError("the fixed points have different Weyl groups"))
  target_length = length(u)
  target_length > length(v) && return zero(t[1])

  R = root_system(W)
  reflections = reflection.(simple_roots(R))
  generator_matrix, _ = _gen_matrix_and_type_of_graph(R)
  T = typeof(one(t[1]))
  terms = Dict{Tuple{typeof(v),Int},T}((one(W), 0) => one(t[1]))
  prefix = one(W)

  for index in word(v)
    i = Int(index)
    # With Oscar's right action this is s₁⋯sⱼ₋₁(αᵢⱼ).
    beta = simple_root(R, i) * inv(prefix)
    beta_class = _billey_root_class(beta, generator_matrix, t)
    next_terms = copy(terms)
    for ((x, k), value) in terms
      k == target_length && continue
      xs = x * reflections[i]
      length(xs) == k + 1 || continue
      key = (xs, k + 1)
      next_terms[key] = get(next_terms, key, zero(value)) + value * beta_class
    end
    terms = next_terms
    prefix *= reflections[i]
  end
  return get(terms, (u, target_length), zero(t[1]))
end

function _billey_root_class(root, generator_matrix, t)
  C = parent(zero(eltype(generator_matrix)))
  coordinates = matrix(C, Oscar.coefficients(root) * generator_matrix)
  return sum(i -> coordinates[i] * t[i], eachindex(t); init=zero(t[1]))
end

function schubert_basis_by_representative(
  ::AbstractGKMGraph;
  representation::Symbol=:localized,
)
  _check_schubert_representation(representation)
  throw(ArgumentError("Billey's formula requires a homogeneous GKM graph with Weyl-group vertices"))
end

function schubert_basis_by_representative(
  G::AbstractGKMGraph{R,V,F};
  representation::Symbol=:localized,
) where {R,V<:GeneralizedFlagVertex,F}
  _check_schubert_representation(representation)
  return Dict(
    label(G, i) => schubert_basis(G, i; representation)
    for i in 1:num_vertices(G)
  )
end
