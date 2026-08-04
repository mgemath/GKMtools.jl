@doc raw"""
    equivariant_coefficient_ring(G)
    equivariant_coefficient_ring(H::GKMCohomology)

Return the polynomial coefficient ring ``H_T^*(pt)`` associated with a GKM
graph or its cohomology parent.
"""
equivariant_coefficient_ring(H::GKMCohomology) = H.coefficient_ring
equivariant_coefficient_ring(G::AbstractGKMGraph) = equivariant_coefficient_ring(get_cohomology(G))

@doc raw"""
    polynomial_gkm_ring(G)

Return the [`Cohomology`](@ref) parent shared by the polynomial and localized
classes on `G`. Polynomial classes are constructed with [`polynomial_class`](@ref).
"""
polynomial_gkm_ring(G::AbstractGKMGraph) = get_cohomology(G)

Base.parent(c::GKMClass) = c.parent
@doc raw"""
    restrictions(c::GKMClass)

Return a copy of the fixed-point restrictions of `c`, ordered by the vertices
of its GKM graph.
"""
Oscar.restrictions(c::GKMClass) = copy(c.restrictions)
Base.length(c::GKMClass) = length(c.restrictions)
Base.getindex(c::GKMClass, v::Integer) = c.restrictions[v]
Base.:(==)(a::GKMClass, b::GKMClass) =
  parent(a) === parent(b) && a.graph === b.graph && a.restrictions == b.restrictions

@doc raw"""
    is_gkm_spline(G, values) -> Bool
    is_gkm_spline(c::GKMClass) -> Bool

Return whether the restrictions satisfy the polynomial GKM relations. For each
edge ``e=(v,w)`` of weight ``α``, this requires
`values[v] - values[w]` to be divisible by ``α`` in ``H_T^*(pt)``.
# Examples
Standard functions for accessing cohomology classes always yield GKM classes:
```jldoctest is_gkm_class
julia> G = projective_space(GKM_graph, 1);

julia> is_gkm_spline(point_class(G, 1))
true
```
Moreover, equivariant cohomology is a ring and a module over the coefficient ring:
```jldoctest is_gkm_class
julia> is_gkm_spline(point_class(G, 1)^2 * point_class(G, 2))
true
```
However, it is possible to cook up non-GKM classes manually.
In the example below, this is because $w(e)=t_1-t_2$, which does not divide $t_1^2 - t_2$.
Here, $e$ is the unique edge of the GKM graph of $\mathbb{P}^1$.
```jldoctest is_gkm_class
julia> (t1, t2) = gens(G.equivariantCohomology.coeffRing);

julia> (e1, e2) = gens(G.equivariantCohomology.cohomRing);

julia> c = t1^2 * e1 + t2 * e2
GKM class with restrictions: 
[t1^2, t2]

julia> is_gkm_spline(c)
false
```
"""
function is_gkm_spline(G::AbstractGKMGraph, values::AbstractVector)
  length(values) == num_vertices(G) || return false
  S = equivariant_coefficient_ring(G)
  local polynomial_values
  try
    polynomial_values = map(values) do value
      try
        S(value)
      catch
        is_unit(denominator(value)) || throw(ArgumentError("non-polynomial restriction"))
        S(divexact(numerator(value), denominator(value)))
      end
    end
  catch
    return false
  end

  t = collect(gens(S))
  for e in edges(G)
    difference = polynomial_values[src(e)] - polynomial_values[dst(e)]
    alpha = _weight_class(G, e, t)
    if iszero(alpha)
      iszero(difference) || return false
    else
      success, _ = divides(difference, alpha)
      success || return false
    end
  end
  return true
end

is_gkm_spline(c::GKMClass) = is_gkm_spline(c.graph, c.restrictions)
# is_gkm_class(args...) = is_gkm_spline(args...)

@doc raw"""
    polynomial_class(G, values; check=true)

Construct a [`GKMClass`](@ref) over ``H_T^*(pt)`` from its fixed-point
restrictions. The length of `values` must equal the number of vertices of `G`.

When `check` is `true`, throw an `ArgumentError` unless the restrictions satisfy
the polynomial GKM edge relations. Set `check=false` only when those relations
are already guaranteed by the caller.
"""
function polynomial_class(G::AbstractGKMGraph, values::AbstractVector; check::Bool=true)
  length(values) == num_vertices(G) || throw(DimensionMismatch(
    "expected $(num_vertices(G)) vertex restrictions, got $(length(values))",
  ))
  polynomial_values = equivariant_coefficient_ring(G).(values)
  check && !is_gkm_spline(G, polynomial_values) && throw(ArgumentError(
    "the restrictions do not satisfy the polynomial GKM edge relations",
  ))
  return GKMClass(get_cohomology(G), G, polynomial_values)
end

@doc raw"""
    localized_class(G, values)

Construct a [`GKMClass`](@ref) whose restrictions lie in the fraction field of
``H_T^*(pt)``. The values are coerced into that field and are not checked for
polynomial GKM divisibility.
"""
function localized_class(G::AbstractGKMGraph, values::AbstractVector)
  length(values) == num_vertices(G) || throw(DimensionMismatch(
    "expected $(num_vertices(G)) vertex restrictions, got $(length(values))",
  ))
  H = get_cohomology(G)
  return GKMClass(H, G, H.localized_coefficient_ring.(values))
end

function _check_same_class_parent(a::GKMClass, b::GKMClass)
  parent(a) === parent(b) && a.graph === b.graph ||
    throw(ArgumentError("GKM classes have different parents"))
  return nothing
end

Base.zero(c::GKMClass) = GKMClass(parent(c), c.graph, zero.(c.restrictions))
Base.one(c::GKMClass) = GKMClass(parent(c), c.graph, one.(c.restrictions))
Base.iszero(c::GKMClass) = all(iszero, c.restrictions)
Base.isone(c::GKMClass) = all(isone, c.restrictions)

function Base.:+(a::GKMClass, b::GKMClass)
  _check_same_class_parent(a, b)
  return GKMClass(parent(a), a.graph, a.restrictions + b.restrictions)
end

function Base.:-(a::GKMClass, b::GKMClass)
  _check_same_class_parent(a, b)
  return GKMClass(parent(a), a.graph, a.restrictions - b.restrictions)
end

Base.:-(a::GKMClass) = GKMClass(parent(a), a.graph, -a.restrictions)

function Base.:*(a::GKMClass, b::GKMClass)
  _check_same_class_parent(a, b)
  return GKMClass(parent(a), a.graph, a.restrictions .* b.restrictions)
end

function Base.:*(f, c::GKMClass)
  scalar = parent(first(c.restrictions))(f)
  return GKMClass(parent(c), c.graph, scalar .* c.restrictions)
end
Base.:*(c::GKMClass, f) = f * c

function Base.:+(c::GKMClass, f)
  scalar = parent(first(c.restrictions))(f)
  return GKMClass(parent(c), c.graph, [x + scalar for x in c.restrictions])
end
Base.:+(f, c::GKMClass) = c + f
Base.:-(c::GKMClass, f) = c + (-f)
Base.:-(f, c::GKMClass) = f + (-c)

function Base.:^(c::GKMClass, n::Integer)
  n >= 0 || throw(DomainError(n, "a GKM class cannot be raised to a negative power"))
  result = one(c)
  factor = c
  exponent = n
  while exponent > 0
    isodd(exponent) && (result *= factor)
    exponent >>= 1
    exponent > 0 && (factor *= factor)
  end
  return result
end

@doc raw"""
    localize(c::GKMClass)

Extend every restriction of `c` to the fraction field of ``H_T^*(pt)``.
If `c` is already localized, return it unchanged.
"""
function localize(c::GKMClass)
  K = parent(c).localized_coefficient_ring
  all(x -> parent(x) === K, c.restrictions) && return c
  return GKMClass(parent(c), c.graph, K.(c.restrictions))
end

function _fraction_to_polynomial(S, f)
  _is_polynomial_fraction(f) || throw(ArgumentError(
    "a fixed-point restriction has a non-polynomial denominator",
  ))
  return S(divexact(numerator(f), denominator(f)))
end

@doc raw"""
    delocalize(G, c::GKMClass; check=true)

Convert a localized class to a polynomial GKM class on `G`.

Throw an `ArgumentError` if a restriction has a genuine denominator. When
`check` is `true`, also verify the polynomial GKM edge relations.
"""
function delocalize(G::AbstractGKMGraph, c::GKMClass; check::Bool=true)
  _check_class_graph(G, c)
  S = equivariant_coefficient_ring(G)
  values = [_fraction_to_polynomial(S, f) for f in c.restrictions]
  return polynomial_class(G, values; check)
end

function _check_class_graph(G::AbstractGKMGraph, c::GKMClass)
  c.graph === G && parent(c) === get_cohomology(G) ||
    throw(ArgumentError("the class belongs to a different GKM graph"))
  return nothing
end

# @doc raw"""
#     is_gkm_class(G, c::GKMClass) -> Bool

# Return whether `c` belongs to the polynomial GKM cohomology of `G`: every
# restriction must be polynomial and all edge-divisibility relations must hold.
# Throw an `ArgumentError` if `c` belongs to another graph.
# """
# function is_gkm_class(G::AbstractGKMGraph, c::GKMClass)
#   _check_class_graph(G, c)
#   S = equivariant_coefficient_ring(G)
#   values = elem_type(S)[]
#   for f in localize(c).restrictions
#     _is_polynomial_fraction(f) || return false
#     push!(values, _fraction_to_polynomial(S, f))
#   end
#   return is_gkm_spline(G, values)
# end

@doc raw"""
    localize_at_vertex(G, c::GKMClass, v::Int)
    localize_at_vertex(G, c::GKMClass, label::String)

Return the fixed-point restriction of `c` at vertex index `v` or at the vertex
with the given `label`. Throw an error if `c` belongs to another graph or the
vertex does not exist.
"""
function localize_at_vertex(G::AbstractGKMGraph, c::GKMClass, v::Int)
  _check_class_graph(G, c)
  checkbounds(c.restrictions, v)
  return c[v]
end

function localize_at_vertex(G::AbstractGKMGraph, c::GKMClass, label::String)
  _check_class_graph(G, c)
  return c[find_vertex_index(label, G)]
end

_localized_vertex_coefficient(G::AbstractGKMGraph, c::GKMClass, v::Int) =
  get_cohomology(G).localized_coefficient_ring(localize_at_vertex(G, c, v))


@doc raw"""
    gens_coeffRing(G)
Return the generators of the equivariant coefficient ring of `G`. The result is a vector of the same length as the torus rank of `G`.
"""
function gens_coeffRing(G)
  return gens(get_cohomology(G).localized_coefficient_ring)
end

@doc raw"""
    gens_cohomRing(G)
Return the generators of the polynomial GKM cohomology ring of `G`. The latter are the classes
dual to the fixed points, with restrictions equal to 1 at one vertex and 0 at all others.
"""
function gens_cohomRing(G::AbstractGKMGraph)
  H = get_cohomology(G)
  K = H.localized_coefficient_ring
  return [GKMClass(H, G, [i == j ? one(K) : zero(K) for i in 1:H.n_vertices])
          for j in 1:H.n_vertices]
end
