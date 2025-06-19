@doc raw"""
    is_gkm_class(c::FreeModElem{QQMPolyRingElem}, G::GKM_graph) -> Bool

Return true if the given class represents an actual cohomology class.
This holds if and only if the difference between localizations at fixed points connected through
an edge e is divisible by the weight of e (see above)

# Examples
Standard functions for accessing cohomology classes always yield GKM classes:
```jldoctest is_gkm_class
julia> G = projective_space(GKM_graph, 1);

julia> is_gkm_class(point_class(G, 1), G)
true
```
Moreover, equivariant cohomology is a ring and a module over the coefficient ring:
```jldoctest is_gkm_class
julia> is_gkm_class(point_class(G, 1)^2 * point_class(G, 2), G)
true
```
However, it is possible to cook up non-GKM classes manually.
In the example below, this is because $w(e)=t_1-t_2$, which does not divide $t_1^2 - t_2$.
Here, $e$ is the unique edge of the GKM graph of $\mathbb{P}^1$.
```jldoctest is_gkm_class
julia> (t1, t2) = gens(G.equivariantCohomology.coeffRing);

julia> (e0, e1) = gens(G.equivariantCohomology.cohomRing);

julia> c = t1^2 * e0 + t2 * e1
t1^2*e[1] + t2*e[2]

julia> is_gkm_class(c, G)
false
```

"""
function is_gkm_class(c::FreeModElem{QQMPolyRingElem}, G::GKM_graph)::Bool

  g = G.g

  for e in edges(g)
    polySrc = c[src(e)]
    polyDst = c[dst(e)]
    dif = polySrc - polyDst
    if !divides(dif, weight_class(e, G))[1]
      return false
    end
  end
  return true
end

function scalar(s::Union{Int, Rational, QQFieldElem}, G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  R = G.equivariantCohomology
  # z = R.coeffRing(s)
  z = G.equivariantCohomology.coeffRing(s)
  return R.cohomRing([z for _ in 1:n_vertices(G.g)])
end

function scalar(s::QQMPolyRingElem, G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  # return R.cohomRing([s for _ in 1:n_vertices(R.gkm.g)])

  return G.equivariantCohomology.cohomRing([s for _ in 1:n_vertices(G.g)])
end 

function zero(G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  return scalar(0,G)
end

function one(G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  return scalar(1,G)
end

@doc raw"""
    point_class(vertexLabel::String, G::GKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the fixed point with given label.

"""
function point_class(vertexLabel::String, G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  @req vertexLabel in G.labels "Vertex not found"
  return point_class(indexin([vertexLabel], G.labels)[1], G)
end

function point_class(G::GKM_graph, vertexLabel::String)::FreeModElem{QQMPolyRingElem}
  return point_class(vertexLabel, G)
end

@doc raw"""
    point_class(vertex::Int, G::GKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual of the given fixed point.

# Examples
```jldoctest point_class
julia> P2 = projective_space(GKM_graph, 2);

julia> point_class(1, P2)
(t1^2 - t1*t2 - t1*t3 + t2*t3)*e[1]

julia> F3 = flag_variety(GKM_graph, [1, 1, 1]);

julia> point_class(1, F3)
(t1^2*t2 - t1^2*t3 - t1*t2^2 + t1*t3^2 + t2^2*t3 - t2*t3^2)*e[1]

```
"""
function point_class(vertex::Int, G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  return point_class(G, vertex)
end 

function point_class(G::GKM_graph, vertex::Int)::FreeModElem{QQMPolyRingElem}
  return euler_class(vertex, G) * gens(G.equivariantCohomology.cohomRing)[vertex]
end 

@doc raw"""
    poincare_dual(gkmSub::GKM_subgraph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant Poincare dual cohomology class of the GKM subgraph.

# Example
```jldoctest poincare_dual
julia> P2 = projective_space(GKM_graph, 2);

julia> P1inP2 = gkm_subgraph_from_vertices(P2, [1, 2])
GKM subgraph of:
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1, 0)

julia> poincare_dual(P1inP2)
(t1 - t3)*e[1] + (t2 - t3)*e[2]
```
"""
function poincare_dual(gkmSub::GKM_subgraph)::FreeModElem{QQMPolyRingElem}

  R = gkmSub.Codomain.equivariantCohomology

  res = zero(gkmSub.Codomain)

  for v in gkmSub.vDict
    vContrib = R.coeffRing(1)
    for w in all_neighbors(gkmSub.Codomain.g, v)
      e = Edge(v,w)
      if !has_edge(gkmSub, e)
        vContrib *= weight_class(e, gkmSub.Codomain)
      end
    end
    res += vContrib * gens(R.cohomRing)[v]
  end

  return res
end

@doc raw"""
    weight_class(e::Edge, G::GKM_graph) -> QQMPolyRingElem

Return the weight of the edge e as an element of the coefficient ring of the equivariant cohomology theory.

# Example
```jldoctest weight_class
julia> H7 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 7))
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (7, 1, 0, -1)
4 -> 1 => (0, 1, 7, -1)
4 -> 3 => (-1, 0, 1, 0)

julia> weight_class(Edge(3, 2), H7)
7*t1 + t2 - t4
```
"""
function weight_class(e::Edge, G::GKM_graph)::QQMPolyRingElem

  R = G.equivariantCohomology
  try
    return R.edgeWeightClasses[e]
  catch err
    if isa(err, KeyError)
      res = _weight_class(e, G)
      R.edgeWeightClasses[e] = res
      return res
    else
      rethrow(err)
    end
  end
end

function _weight_class(e::Edge, G::GKM_graph)::QQMPolyRingElem
  R = G.equivariantCohomology
  w = G.w[e]
  rk = rank_torus(G)
  coeffs = R.coeffRing
  t = gens(coeffs)
  
  res = coeffs(0)
  for i in 1:rk
    res += w[i] * t[i]
  end
  return res
end

@doc raw"""
    euler_class(vertex::Int64, G::GKM_graph) -> QQMPolyRingElem

Return the Euler class of the normal bundle of the fixed point v.
This is the localization of `point_class(vertex, G)` to the given vertex.

# Example
```jldoctest euler_class
julia> F3 = flag_variety(GKM_graph, [1, 1, 1]);

julia> euler_class(1, F3)
t1^2*t2 - t1^2*t3 - t1*t2^2 + t1*t3^2 + t2^2*t3 - t2*t3^2

julia> point_class(1, F3)
(t1^2*t2 - t1^2*t3 - t1*t2^2 + t1*t3^2 + t2^2*t3 - t2*t3^2)*e[1]
```
"""
function euler_class(vertex::Int64, G::GKM_graph)::QQMPolyRingElem

  R = G.equivariantCohomology
  res = R.pointEulerClasses[vertex]
  if isnothing(res)
    res = _euler_class(vertex, G)
    R.pointEulerClasses[vertex] = res
  end
  return res
end

function _euler_class(vertex::Int,  G::GKM_graph)::QQMPolyRingElem

  R = G.equivariantCohomology
  res = R.coeffRing(1)
  for i in all_neighbors(G.g, vertex)
    res = mul!(res, weight_class(Edge(vertex, i), G))
  end
  return res
end

@doc raw"""
    euler_class(vertexLabel::String, G::GKM_graph) -> QQMPolyRingElem

Return the euler class of the normal bundle of the fixed point with the given label.
"""
function euler_class(vertexLabel::String, G::GKM_graph)::QQMPolyRingElem
  @req vertexLabel in G.labels "Vertex not found"
  return euler_class(indexin([vertexLabel], G.labels)[1], G)
end

@doc raw"""
    integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, G::GKM_graph; check::Bool=true) -> QQMPolyRingElem

Integrate the GKM class, yielding an element of the coefficient ring. This checks if `is_gkm_class(class,R) == true` and throws an error otherwise.

# Examples
```jldoctest integrate_gkm_class
julia> P2 = projective_space(GKM_graph, 2);

julia> integrate_gkm_class(point_class(1, P2), P2)
1

julia> P2inP1 = gkm_subgraph_from_vertices(P2, [1, 2]);

julia> pd = poincare_dual(P2inP1);

julia> integrate_gkm_class(pd, P2)
0

julia> integrate_gkm_class(pd^2, P2)
1

julia> (t1, t2, t3) = gens(P2.equivariantCohomology.coeffRing);

julia> integrate_gkm_class(t3 * pd^2 + (t2^2 - t1)*point_class(3, P2), P2)
-t1 + t2^2 + t3
```
"""
function integrate_gkm_class(class::FreeModElem{QQMPolyRingElem}, G::GKM_graph; check::Bool=true)::QQMPolyRingElem
  R = G.equivariantCohomology

  if check
    @req is_gkm_class(class, G) "Given class is not GKM."
  end
  res = integrate(class, G)
  (flag, quotient) = divides(numerator(res), denominator(res))
  if check
    @req flag "Integral of GKMclass is fraction"
  end
  return quotient
end


@doc raw"""
    first_chern_class(G::GKM_graph) -> FreeModElem{QQMPolyRingElem}

Return the equivariant first Chern class of the GKM space, i.e., $c_1^T(TX)$ where $TX$ is the ($T$-equivariant) tangent bundle of $X$.

# Examples
```jldoctest first_chern_class
julia> P2 = projective_space(GKM_graph, 2);

julia> first_chern_class(P2)
(2*t1 - t2 - t3)*e[1] + (-t1 + 2*t2 - t3)*e[2] + (-t1 - t2 + 2*t3)*e[3]

julia> H3 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 3))
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (3, 1, 0, -1)
4 -> 1 => (0, 1, 3, -1)
4 -> 3 => (-1, 0, 1, 0)

julia> first_chern_class(H3)
(-t1 - t2 - 2*t3 + t4)*e[1] + (-2*t1 - t2 - t3 + t4)*e[2] + (4*t1 + t2 - t3 - t4)*e[3] + (-t1 + t2 + 4*t3 - t4)*e[4]
```
"""
function first_chern_class(G::GKM_graph)::FreeModElem{QQMPolyRingElem}

  R = G.equivariantCohomology
  res = zero(G)
  for v in 1:n_vertices(G.g)
    localFactor = zero(R.coeffRing)
    for w in all_neighbors(G.g, v)
      localFactor += weight_class(Edge(v, w), G)
    end
    res += localFactor * gens(R.cohomRing)[v]
  end
  return res
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, G::GKM_graph, e::Edge) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class over the curve represented by the GKM graph edge.
In mathematical notation, if the edge is $e = (p\rightarrow q)$ and the class is $c=(f_v)_{v\in X^T}$ then
$\int_e c = \frac{f_p - f_q}{w(e)}$ where $w(e)$ is the weight of $e$.
If `is_gkm_class(class, G) == true` then the result is a polynomial in the variables `gens(G.equivariantCohomology.coeffRing)`.

# Example
```jldoctest integrate_edge
julia> P2 = projective_space(GKM_graph, 2);

julia> integrate(first_chern_class(P2), P2, Edge(1, 2))
3
```
"""
function integrate(
  class::FreeModElem{QQMPolyRingElem},
  G::GKM_graph,
  e::Edge
)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  return (class[src(e)] - class[dst(e)]) // weight_class(e, G)
end

@doc raw"""
    integrate(class::FreeModElem{QQMPolyRingElem}, G::GKM_graph) -> AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

Integrate the cohomology class, yielding an element of the fraction field of the equivariant coefficient ring.

This uses the Atiyah--Bott localization formula:  If $c = (f_v)_{v\in X^T}$ then $\int_X c = \sum_{v\in X^T}\frac{f_v}{e_T(T_vX)}$
where $e_T(T_vX)$ is the equivariant Euler class of the tangent bundle of $X$ at the fixed point $v$.

If `is_gkm_class(class, G) == true`, this fraction will be a polynomial.
!!! note
    Use `integrate_gkm_class()` instead if you know that the class is a GKM class.

# Examples
```jldoctest integrate_global
julia> G24 = flag_variety(GKM_graph, [2, 4]); # Grassmannian of 2-planes in C^4

julia> integrate(point_class(G24, 1), G24)
1
```
In contrast to `integrate_gkm_class`, we can also integrate tuples $(f_v)_{v\in X^T}$ that do not satisfy `is_gkm_class(class) == true`:
```jldoctest integrate_global
julia> P1 = projective_space(GKM_graph, 1);

julia> (t1, t2) = gens(P1.equivariantCohomology.coeffRing);

julia> (e0, e1) = gens(P1.equivariantCohomology.cohomRing);

julia> c = t1^2 * e0 + t2 * e1
t1^2*e[1] + t2*e[2]

julia> is_gkm_class(c, P1)
false

julia> integrate(c, P1)
(t1^2 - t2)//(t1 - t2)
```
"""
function integrate(class::FreeModElem{QQMPolyRingElem}, G::GKM_graph)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  
  R = G.equivariantCohomology
  res = R.coeffRing(0)
  for v in 1:n_vertices(G.g)
    euler = euler_class(v, G)
    contrib = class[v] // euler
    res = add!(res, contrib)
  end
  return res
end

function multiply(a::FreeModElem{QQMPolyRingElem}, b::FreeModElem{QQMPolyRingElem}, G::GKM_graph)::FreeModElem{QQMPolyRingElem}
  return R.cohomRing([a[i] * b[i] for i in 1:n_vertices(G.g)])
end

function *(a::FreeModElem{T}, b::FreeModElem{T})::FreeModElem{T} where T <: RingElem
  return parent(a)([a[i] * b[i] for i in 1:rank(parent(a))])
end

function ^(a::FreeModElem{T}, n::Int64)::FreeModElem{T} where T <: RingElem
  return parent(a)([a[i]^n for i in 1:rank(parent(a))])
end

function Base.show(io::IO, R::GKM_cohomology_ring)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM cohomology ring")
  else
    # nested printing allowed, preferably terse
    print(io, "Cohomology ring of GKM graph with $(n_vertices(R.gkm.g)) nodes, valency $(valency(R.gkm)), and rank $(rank_torus(R.gkm)) torus")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", R::GKM_cohomology_ring)

  print(io, "Cohomology ring of GKM graph with $(n_vertices(R.gkm.g)) nodes, valency $(valency(R.gkm)), and rank $(rank_torus(R.gkm)) torus\n")
  print(io, "Coefficient Ring: ")
  show(io, MIME"text/plain"(), R.coeffRing)
  print(io, "\nCohomology Ring: ")
  show(io, MIME"text/plain"(), R.cohomRing)
  
end