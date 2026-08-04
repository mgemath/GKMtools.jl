# Cohomology

Let $X$ be a GKM space with respect to the complex torus $T$.
By [GKM98](@cite) and [GZ01; Theorem 1.7.3](@cite), we have the following description of its equivariant cohomology ring $H^*_T(X;\mathbb{Q})$.

Each fixed point $x\in X$ gives a ring homomorphism $H^*_T(X;\mathbb{Q})\rightarrow H^*_T(\{x\};\mathbb{Q})\cong \mathbb{Q}[\mathfrak{t}]$.
We may combine these maps by taking all fixed points at once.
Let $V$ be the set of fixed points of $X$ and $E$ the set of $T$-invariant rational curves.
That is, $V$ are the vertices of the GKM graph and $E$ the edges.
Then the map

```math
    H^*_T(X;\mathbb{Q}) \longrightarrow \bigoplus_{x\in V} \mathbb{Q}[\mathfrak{t}]
```

is injective and its image consists of all $(f_x)_{x\in V}$ such that $f_{\text{src}(e)}\equiv f_{\text{dst}(e)}$ mod $w(e)$ for all edges $e\in E$, where $w(e)$ is the weight of $e$.

We may further identify $\mathbb{Q}[\mathfrak{t}]\cong \mathbb{Q}[t_1,\dots,t_r]$ where $r=\dim_{\mathbb{C}}(T)$.
Hence, this package represents elements of $H^*_T(X;\mathbb{Q})$ as tuples of polynomials indexed by vertices of the GKM graph.

The function `equivariant_coefficient_ring` returns the coefficient ring of the cohomology ring.

```jldoctest cohomology
julia> G = projective_space(GKMGraph, 2);

julia> S = equivariant_coefficient_ring(G) # == \mathbb{Q}[t_1,\dots,t_r] == H^*_T(\{x\};\mathbb{Q})
Multivariate polynomial ring in 3 variables t1, t2, t3
  over rational field
```
In order to construct an element of the cohomology ring of `G`, we first define a vector of elements of `S` of length 
`num_vertices(G)`, such that each entry is the localization of that class to the vertex number. 
After that, we apply `polynomial_class` to that vector. This will check if the vector represents an 
element of the cohomology ring. If this check fails, we get an error. For example, the elements zero and one of the cohomology ring can be show as follows:

```jldoctest cohomology
julia> v0 = [zero(S) for _ in 1:num_vertices(G)] # we constructed the vector
3-element Vector{QQMPolyRingElem}:
 0
 0
 0

julia> polynomial_class(G, v0)
GKM class with restrictions: 
[0, 0, 0]

julia> v1 = [one(S) for _ in 1:num_vertices(G)]
3-element Vector{QQMPolyRingElem}:
 1
 1
 1

julia> polynomial_class(G, v1)
GKM class with restrictions: 
[1, 1, 1]
```
Taking a vector that does not satisfy the condition  $f_{\text{src}(e)}\equiv f_{\text{dst}(e)}$ mod $w(e)$ leads to an error.

```julia-repl
julia> v2 = [one(S) for _ in 1:num_vertices(G)];

julia> v2[1] = zero(S);

julia> polynomial_class(G, v2)
ERROR: ArgumentError: the restrictions do not satisfy the polynomial GKM edge relations
Stacktrace:
```


## General functions
```@docs
GKMClass
is_gkm_spline
Oscar.restrictions
weight_class
gens_coeffRing
gens_cohomRing
unit_cohomology_ring
localize_at_vertex
localize
delocalize
integrate
integrate_gkm_class
equivariant_coefficient_ring
polynomial_gkm_ring
polynomial_class
```

## Subvarieties
```@docs
point_class
poincare_dual
```

## Chern Classes
```@docs
first_chern_class
chern_class
GKMtools.total_chern_class
chern_classes
```