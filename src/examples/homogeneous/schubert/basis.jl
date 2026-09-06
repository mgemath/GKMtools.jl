@doc raw"""
    schubert_basis(G::AbstractGKMGraph; representation=:localized) -> Vector

Return Billey's equivariant Schubert basis of a homogeneous GKM variety
`G = G/P`. The result is grouped by complex codimension: `result[d + 1]`
contains the classes of codimension `d` (cohomological degree `2d`). Within a
group, classes follow the fixed-point ordering of `G`.

Set `representation=:polynomial` to return [`GKMClass`](@ref) objects over
`H_T^*(pt)`. The default `:localized` preserves the historical fraction-field
representation.
"""
function schubert_basis(
  G::AbstractGKMGraph{R,V,F};
  representation::Symbol=:localized,
) where {R,V<:GeneralizedFlagVertex,F}
  _check_schubert_representation(representation)
  vertex_flags = flag.(vertices_structure(G))
  dimensions = length.(vertex_flags)
  t = collect(gens(equivariant_coefficient_ring(G)))
  zero_value = zero(first(t))

  # Compute all Schubert restrictions at a fixed point in one traversal of
  # its reduced word. The previous implementation repeated this traversal
  # once for every Schubert class, making the whole basis quadratic in the
  # number of fixed points before accounting for the subword recursion.
  values = [fill(zero_value, num_vertices(G)) for _ in vertex_flags]
  flag_to_index = Dict(u => i for (i, u) in enumerate(vertex_flags))
  allowed_flags = _right_descent_closure(vertex_flags)
  for (fixed_index, v) in enumerate(vertex_flags)
    for (u, value) in _billey_localizations(G, v, t, allowed_flags)
      target_index = get(flag_to_index, u, 0)
      iszero(target_index) || (values[target_index][fixed_index] = value)
    end
  end

  polynomial_classes = [polynomial_class(G, row; check=false) for row in values]
  classes = representation == :polynomial ? polynomial_classes : localize.(polynomial_classes)
  basis = [eltype(classes)[] for _ in 0:maximum(dimensions)]
  for vertex in eachindex(classes)
    push!(basis[dimensions[vertex] + 1], classes[vertex])
  end
  return basis
end

function _right_descent_closure(elements)
  isempty(elements) && return Set(elements)
  R = root_system(parent(first(elements)))
  reflections = reflection.(simple_roots(R))
  closure = Set(elements)
  queue = collect(elements)
  head = 1
  while head <= length(queue)
    element = queue[head]
    head += 1
    for reflection in reflections
      predecessor = element * reflection
      length(predecessor) + 1 == length(element) || continue
      predecessor in closure && continue
      push!(closure, predecessor)
      push!(queue, predecessor)
    end
  end
  return closure
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
  1 <= vertex <= num_vertices(G) || throw(BoundsError(vertices_structure(G), vertex))

  u = flag(vertices_structure(G)[vertex])
  t = collect(gens(equivariant_coefficient_ring(G)))
  root_data = root_system(parent(u))
  reflections = reflection.(simple_roots(root_data))
  lower_interval = _bruhat_lower_interval(u, reflections)
  values = [
    _billey_localization(
      G, u, flag(vertices_structure(G)[v]), t, lower_interval,
    )
    for v in 1:num_vertices(G)
  ]
  # Billey's formula produces a GKM spline by construction. Rechecking all
  # edge divisibilities is especially expensive for large homogeneous graphs.
  polynomial = polynomial_class(G, values; check=false)
  return representation == :polynomial ? polynomial : localize(polynomial)
end

function _bruhat_lower_interval(u, reflections)
  W = parent(u)
  interval = Set([one(W)])
  for index in word(u)
    reflection = reflections[Int(index)]
    previous = collect(interval)
    for x in previous
      xs = x * reflection
      length(xs) == length(x) + 1 || continue
      push!(interval, xs)
    end
  end
  return interval
end

# Return Billey's localizations for every reduced subword of v. Since the
# minimal parabolic representatives form a lower Bruhat ideal, at most one
# state per vertex of G survives, even for a large ambient Weyl group.
function _billey_localizations(
  G::AbstractGKMGraph,
  v,
  t=gens_coeffRing(G),
  allowed=nothing,
)
  W = parent(v)
  R = root_system(W)
  reflections = reflection.(simple_roots(R))
  generator_matrix, _ = _gen_matrix_and_type_of_graph(R)
  T = typeof(one(t[1]))
  terms = Dict{typeof(v),T}(one(W) => one(t[1]))
  prefix = one(W)

  for index in word(v)
    i = Int(index)
    beta = simple_root(R, i) * inv(prefix)
    beta_class = _billey_root_class(beta, generator_matrix, t)
    next_terms = copy(terms)
    for (x, value) in terms
      xs = x * reflections[i]
      length(xs) == length(x) + 1 || continue
      isnothing(allowed) || xs in allowed || continue
      next_terms[xs] = get(next_terms, xs, zero(value)) + value * beta_class
    end
    terms = next_terms
    prefix *= reflections[i]
  end
  return terms
end

# Dynamic programming collects the root products of reduced subwords, avoiding
# an explicit enumeration of all 2^length(v) subwords.
function _billey_localization(
  G::AbstractGKMGraph,
  u,
  v,
  t=gens_coeffRing(G),
  lower_interval=nothing,
)
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

      # A selected prefix can contribute to u only when it lies below u in
      # Bruhat order. Without this test we retain nearly every short subword
      # of v, which is prohibitive for exceptional Weyl groups.
      (isnothing(lower_interval) ? (xs == u || xs < u) : (xs in lower_interval)) || continue
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

@doc raw"""
    schubert_class(schubert::AbstractGKMSubgraph, label::String) -> GKMClass

Return the equivariant Poincaré dual, inside `schubert`, of the Schubert
subvariety indexed by `label`. The lower interval is determined directly from
the Weyl-group element stored in the `flag` field of each vertex; no separately
constructed Bruhat order is required.

# Example
```jldoctest schubert_classes
julia> R = root_system(:A, 2);

julia> S = generalized_gkm_schubert(R, "s1*s2");

julia> schubert_class(S, "s1*s2")
GKM class with restrictions: 
[1, 1, 1, 1]

julia> schubert_class(S, "s1")
GKM class with restrictions: 
[t2 - t3, t1 - t3, 0, 0]
```
"""
function schubert_class(
  schubert::AbstractGKMSubgraph,
  vertex_label::String,
)
  G = subgraph(schubert)
  target_index = findfirst(v -> label(G, v) == vertex_label, vertices(G))
  isnothing(target_index) && throw(ArgumentError(
    "Schubert vertex label not found: $vertex_label",
  ))

  target_flag = vertices_structure(G)[target_index].flag
  interval = [
    v for v in vertices(G)
    if let vertex_flag = vertices_structure(G)[v].flag
      vertex_flag == target_flag || vertex_flag < target_flag
    end
  ]
  return poincare_dual(subgraph_from_vertices(G, interval))
end

@doc raw"""
    schubert_classes(schubert::AbstractGKMSubgraph) -> Matrix

Return the fixed-point restrictions of all Schubert classes on `schubert`.
Row `v` contains the restrictions of the class indexed by the `v`-th vertex.

# Example
```jldoctest schubert_classes
julia> schubert_classes(S)
4×4 Matrix{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}:
 t1*t2 - t1*t3 - t2^2 + t2*t3  0        0        0
 t2 - t3                       t1 - t3  0        0
 t1 - t2                       0        t1 - t2  0
 1                             1        1        1
```
"""
function Oscar.schubert_classes(schubert::AbstractGKMSubgraph)
  G = subgraph(schubert)
  classes = [schubert_class(schubert, label(G, v)) for v in vertices(G)]
  values = Matrix{eltype(restrictions(first(classes)))}(
    undef, num_vertices(G), num_vertices(G),
  )
  for v in vertices(G)
    values[v, :] = restrictions(classes[v])
  end
  return values
end
