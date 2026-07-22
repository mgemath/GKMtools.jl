"""
    GKMSplineRing(G)

The polynomial GKM (or spline) ring of `G` over
`S = H_T^*(pt)`. Its elements are polynomial tuples `(f_v)` satisfying
`f_v - f_w ∈ (α_vw)` on every edge. No freeness, basis, or geometric
realizability assumption is made.
"""
struct GKMSplineRing{G,S}
  graph::G
  coefficient_ring::S
end

"""An element of a [`GKMSplineRing`](@ref), stored by vertex restrictions."""
struct GKMSpline{P,R}
  parent::P
  restrictions::Vector{R}
end

equivariant_coefficient_ring(H::GKMCohomology) = H.coefficient_ring
equivariant_coefficient_ring(G::AbstractGKMGraph) = equivariant_coefficient_ring(get_cohomology(G))

polynomial_gkm_ring(G::AbstractGKMGraph) =
  GKMSplineRing(G, equivariant_coefficient_ring(G))

Base.parent(c::GKMSpline) = c.parent
AbstractAlgebra.coefficient_ring(H::GKMSplineRing) = H.coefficient_ring
graph(H::GKMSplineRing) = H.graph
Oscar.restrictions(c::GKMSpline) = copy(c.restrictions)
Base.length(c::GKMSpline) = length(c.restrictions)
Base.getindex(c::GKMSpline, v::Integer) = c.restrictions[v]

Base.:(==)(H::GKMSplineRing, K::GKMSplineRing) =
  H.graph === K.graph && H.coefficient_ring === K.coefficient_ring
Base.:(==)(a::GKMSpline, b::GKMSpline) =
  parent(a) == parent(b) && a.restrictions == b.restrictions

function Base.show(io::IO, H::GKMSplineRing)
  print(io, "polynomial GKM spline ring on $(num_vertices(H.graph)) vertices over ")
  show(io, H.coefficient_ring)
end

function Base.show(io::IO, c::GKMSpline)
  print(io, "GKM spline with restrictions ")
  show(io, c.restrictions)
end

function _spline_edge_weight(H::GKMSplineRing, e::Edge)
  return _weight_class(H.graph, e, collect(gens(H.coefficient_ring)))
end

"""
    is_gkm_spline(G, restrictions) -> Bool
    is_gkm_spline(c::GKMSpline) -> Bool

Test the polynomial edge-divisibility relations. This does not assert that the
spline ring is the equivariant cohomology of a geometric space.
"""
function is_gkm_spline(G::AbstractGKMGraph, values::AbstractVector)
  length(values) == num_vertices(G) || return false
  H = polynomial_gkm_ring(G)
  S = coefficient_ring(H)
  local polynomial_values
  try
    polynomial_values = S.(values)
  catch
    return false
  end

  for e in edges(G)
    difference = polynomial_values[src(e)] - polynomial_values[dst(e)]
    alpha = _spline_edge_weight(H, e)
    if iszero(alpha)
      iszero(difference) || return false
    else
      success, _ = divides(difference, alpha)
      success || return false
    end
  end
  return true
end

is_gkm_spline(c::GKMSpline) = is_gkm_spline(graph(parent(c)), c.restrictions)

"""
    polynomial_class(H, restrictions; check=true)
    polynomial_class(G, restrictions; check=true)

Construct a polynomial GKM spline. Public construction checks all edge
relations by default; internal callers known to produce splines may set
`check=false`.
"""
function polynomial_class(H::GKMSplineRing, values::AbstractVector; check::Bool=true)
  length(values) == num_vertices(graph(H)) || throw(DimensionMismatch(
    "expected $(num_vertices(graph(H))) vertex restrictions, got $(length(values))",
  ))
  polynomial_values = coefficient_ring(H).(values)
  check && !is_gkm_spline(graph(H), polynomial_values) && throw(ArgumentError(
    "the restrictions do not satisfy the polynomial GKM edge relations",
  ))
  return GKMSpline(H, polynomial_values)
end

polynomial_class(G::AbstractGKMGraph, values::AbstractVector; check::Bool=true) =
  polynomial_class(polynomial_gkm_ring(G), values; check)

function _check_same_spline_parent(a::GKMSpline, b::GKMSpline)
  parent(a) == parent(b) || throw(ArgumentError("GKM splines have different parents"))
  return nothing
end

Base.zero(H::GKMSplineRing) =
  GKMSpline(H, fill(zero(coefficient_ring(H)), num_vertices(graph(H))))
Base.one(H::GKMSplineRing) =
  GKMSpline(H, fill(one(coefficient_ring(H)), num_vertices(graph(H))))
Base.zero(c::GKMSpline) = zero(parent(c))
Base.one(c::GKMSpline) = one(parent(c))
Base.iszero(c::GKMSpline) = all(iszero, c.restrictions)
Base.isone(c::GKMSpline) = all(isone, c.restrictions)

function Base.:+(a::GKMSpline, b::GKMSpline)
  _check_same_spline_parent(a, b)
  return GKMSpline(parent(a), a.restrictions + b.restrictions)
end

function Base.:-(a::GKMSpline, b::GKMSpline)
  _check_same_spline_parent(a, b)
  return GKMSpline(parent(a), a.restrictions - b.restrictions)
end

Base.:-(a::GKMSpline) = GKMSpline(parent(a), -a.restrictions)

function Base.:*(a::GKMSpline, b::GKMSpline)
  _check_same_spline_parent(a, b)
  return GKMSpline(parent(a), a.restrictions .* b.restrictions)
end

function Base.:*(f, c::GKMSpline)
  scalar = coefficient_ring(parent(c))(f)
  return GKMSpline(parent(c), scalar .* c.restrictions)
end
Base.:*(c::GKMSpline, f) = f * c

function Base.:+(c::GKMSpline, f)
  scalar = coefficient_ring(parent(c))(f)
  return GKMSpline(parent(c), [x + scalar for x in c.restrictions])
end
Base.:+(f, c::GKMSpline) = c + f
Base.:-(c::GKMSpline, f) = c + (-f)
Base.:-(f, c::GKMSpline) = f + (-c)

function Base.:^(c::GKMSpline, n::Integer)
  n >= 0 || throw(DomainError(n, "a polynomial GKM spline cannot be raised to a negative power"))
  result = one(parent(c))
  factor = c
  exponent = n
  while exponent > 0
    isodd(exponent) && (result *= factor)
    exponent >>= 1
    exponent > 0 && (factor *= factor)
  end
  return result
end

"""Embed a polynomial GKM spline in the existing fraction-field localization."""
function localize(c::GKMSpline)
  H = parent(c)
  G = graph(H)
  localized = get_cohomology(G).localized_cohomology
  K = get_cohomology(G).localized_coefficient_ring
  e = gens_cohomRing(G)
  result = zero(localized)
  for v in vertices(G)
    iszero(c[v]) || (result += K(c[v]) * e[v])
  end
  return result
end

function _fraction_to_polynomial(S, f)
  _is_polynomial_fraction(f) || throw(ArgumentError(
    "a fixed-point restriction has a non-polynomial denominator",
  ))
  return divexact(numerator(f), denominator(f))
end

"""
    delocalize(G, c; check=true)

Convert a localized class to a polynomial spline. This fails when a restriction
has a genuine denominator or the polynomial restrictions violate an edge
relation.
"""
function delocalize(G::AbstractGKMGraph, c; check::Bool=true)
  localized = get_cohomology(G).localized_cohomology
  parent(c) == localized || throw(ArgumentError(
    "the class does not belong to the localized cohomology ring of G",
  ))
  S = equivariant_coefficient_ring(G)
  values = [
    S(_fraction_to_polynomial(S, _localized_vertex_coefficient(G, c, v)))
    for v in vertices(G)
  ]
  return polynomial_class(G, values; check)
end

function _check_spline_graph(G::AbstractGKMGraph, c::GKMSpline)
  graph(parent(c)) === G || throw(ArgumentError("the spline belongs to a different GKM graph"))
  return nothing
end

function is_gkm_class(G::AbstractGKMGraph, c::GKMSpline)
  _check_spline_graph(G, c)
  return is_gkm_spline(c)
end

function localize_at_vertex(G::AbstractGKMGraph, c::GKMSpline, v::Int)
  _check_spline_graph(G, c)
  checkbounds(c.restrictions, v)
  return c[v]
end

function localize_at_vertex(G::AbstractGKMGraph, c::GKMSpline, label::String)
  _check_spline_graph(G, c)
  return c[find_vertex_index(label, G)]
end

function integrate(G::AbstractGKMGraph, c::GKMSpline, e::Edge)
  _check_spline_graph(G, c)
  return integrate(G, localize(c), e)
end
function integrate(G::AbstractGKMGraph, c::GKMSpline, s1::String, s2::String)
  _check_spline_graph(G, c)
  return integrate(G, localize(c), s1, s2)
end
function integrate(G::AbstractGKMGraph, c::GKMSpline)
  _check_spline_graph(G, c)
  return integrate(G, localize(c))
end
