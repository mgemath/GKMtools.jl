
function _euler_class(G, v::Int, t::Vector{T}) where {T}
  res = one(t[1])
  for i in eachindex(flags(G, v))
    res *= _flag_weight_class(G, v, i, t)
  end
  return res
end

function _euler_class(G, v::Int)
  return _euler_class(G, v, gens_coeffRing(G))
end


@doc raw"""
    point_class(G::AbstractGKMGraph, v::Int)
    point_class(v::Int, G::AbstractGKMGraph)
    point_class(G::AbstractGKMGraph, label::String)
    point_class(label::String, G::AbstractGKMGraph)

Return the equivariant cohomology class Poincaré dual to the fixed point at
vertex `v`, or at the vertex with the given `label`. Its only nonzero
restriction is the equivariant Euler class of the tangent representation at
that vertex.

Throw an error when the vertex index or label does not exist.
"""
point_class(
  v::Int,
  G::AbstractGKMGraph,
) = point_class(G, v)

function point_class(G::AbstractGKMGraph, v::Int)
  checkbounds(vertices_structure(G), v)
  return _euler_class(G, v) * gens_cohomRing(G)[v]
end

point_class(
  vertex_label::String,
  G::AbstractGKMGraph,
) = point_class(G, vertex_label)

function point_class(
  G::AbstractGKMGraph,
  vertex_label::String,
)
  index = find_vertex_index(vertex_label, G)
  return point_class(G, index)
end

_is_polynomial_fraction(f) = is_unit(denominator(f))


@doc raw"""
    integrate(G::AbstractGKMGraph, c)
    integrate(G::AbstractGKMGraph, c, e::Edge)
    integrate(G::AbstractGKMGraph, c, source::String, destination::String)

Integrate the GKM class `c` by equivariant localization.

With no edge argument, return the Atiyah--Bott sum over all vertices,
```math
    \sum_{v\in V} \frac{c|_v}{e_T(T_vG)}.
```
With an edge or its endpoint labels, integrate over
the corresponding invariant curve and return the quotient of `(c|_v-c|_w)` and 
`weight_class(G,e)`.

The result belongs to the localized equivariant coefficient ring. Throw an
`ArgumentError` if the supplied edge is not an edge of `G`.
"""
function integrate(
  G::AbstractGKMGraph,
  c,
  e::Edge,
)
  has_edge(graph(G), src(e), dst(e)) ||
    throw(ArgumentError("$e is not an edge of G"))

  c_v = _localized_vertex_coefficient(G, c, src(e))
  c_w = _localized_vertex_coefficient(G, c, dst(e))

  return (c_v - c_w) * inv(weight_class(G, e))
end

function integrate(G, c, s1::String, s2::String)
  e = Edge(find_vertex_index(s1, G), find_vertex_index(s2, G))
  return integrate(G, c, e)
end

function integrate(G, c)
  t = gens_coeffRing(G)
  answer = zero(parent(t[1]))

  for v in vertices(G)
    localization = _localized_vertex_coefficient(G, c, v)
    euler = _euler_class(G, v, t)
    answer += localization * inv(euler)
  end

  return answer
end
