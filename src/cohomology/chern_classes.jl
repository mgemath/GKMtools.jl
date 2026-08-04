
function _euler_class(G::AbstractGKMGraph, v::Int, t::Vector{T}) where {T}
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

# Examples
```jldoctest point_class
julia> P2 = projective_space(GKMGraph, 2);

julia> point_class(1, P2)
GKM class with restrictions: 
[t1^2 - t1*t2 - t1*t3 + t2*t3, 0, 0]

julia> F3 = flag_variety(GKMGraph, [1, 1, 1]);

julia> point_class(1, F3)
GKM class with restrictions: 
[t1^2*t2 - t1^2*t3 - t1*t2^2 + t1*t3^2 + t2^2*t3 - t2*t3^2, 0, 0, 0, 0, 0]
```
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
    integrate(c)
    integrate(c, e::Edge)
    integrate(c, source::String, destination::String)

Integrate the GKM class `c` by equivariant localization.

With no edge argument, return the Atiyah--Bott sum over all vertices,
```math
    \sum_{v\in V} \frac{c|_v}{e_T(T_vG)}.
```
With an edge or its endpoint labels, integrate over
the corresponding invariant curve and return the quotient of `(c|_v-c|_w)` and 
`weight_class(G,e)`.

The result belongs to the localized equivariant coefficient ring. Throw an
`ArgumentError` if the supplied edge is not an edge of `graph(c)`.

# Example
```jldoctest integrate_edge
julia> P2 = projective_space(GKMGraph, 2);

julia> integrate(first_chern_class(P2), Edge(1, 2))
3
```
In contrast to `integrate_gkm_class`, we can also integrate tuples $(f_v)_{v\in X^T}$ that do not satisfy `is_gkm_spline(c)`. For example, the following class is not a GKM spline, but we can still integrate it over the edge from vertex 1 to vertex 2:
```jldoctest integrate_global
julia> P1 = projective_space(GKMGraph, 1);

julia> (t1, t2) = gens(P1.equivariantCohomology.coeffRing);

julia> (e1, e2) = gens(P1.equivariantCohomology.cohomRing);

julia> c = t1^2 * e1 + t2 * e2
GKM class with restrictions: 
[t1^2, t2]

julia> is_gkm_spline(c)
false

julia> integrate(c)
(t1^2 - t2)//(t1 - t2)
```
"""
function integrate(
  c,
  e::Edge,
)
  G = graph(c)
  has_edge(graph(G), src(e), dst(e)) ||
    throw(ArgumentError("$e is not an edge of G"))

  c_v = _localized_vertex_coefficient(G, c, src(e))
  c_w = _localized_vertex_coefficient(G, c, dst(e))

  return (c_v - c_w) * inv(weight_class(G, e))
end

function integrate(c, s1::String, s2::String)
  G = graph(c)
  e = Edge(find_vertex_index(s1, G), find_vertex_index(s2, G))
  return integrate(c, e)
end

function integrate(c)
  G = graph(c)
  t = gens_coeffRing(G)
  answer = zero(parent(t[1]))

  for v in vertices(G)
    localization = _localized_vertex_coefficient(G, c, v)
    euler = _euler_class(G, v, t)
    answer += localization * inv(euler)
  end

  return answer
end

@doc raw"""
    integrate_gkm_class(c)

Integrate the GKM class, yielding an element of the coefficient ring. This checks if `is_gkm_spline(c)` is true and throws an error otherwise.

# Examples
```jldoctest integrate_gkm_class
julia> P2 = projective_space(GKMGraph, 2);

julia> integrate_gkm_class(point_class(1, P2))
1

julia> P2inP1 = subgraph_from_vertices(P2, [1, 2]);

julia> pd = poincare_dual(P2inP1);

julia> integrate_gkm_class(pd)
0

julia> integrate_gkm_class(pd^2)
1

julia> (t1, t2, t3) = gens_coeffRing(P2);

julia> integrate_gkm_class(t3 * pd^2 + (t2^2 - t1)*point_class(3, P2))
-t1 + t2^2 + t3
```
"""
function integrate_gkm_class(c::GKMClass)
  is_gkm_spline(c) || throw(ArgumentError("the class is not a GKM spline"))
  return integrate(c)
end