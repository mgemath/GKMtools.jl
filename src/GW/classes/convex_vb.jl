struct _ConvexBundleFunctor end
struct _ConcaveBundleFunctor end

@inline _decorated_graph(dt::GW_decorated_tree) = dt.tree
@inline _decorated_graph(dt::GW_decorated_graph) = dt.g
@inline _edge_degrees(dt::GW_decorated_tree, e::Edge) = (edgeMult(e, dt),)
@inline _edge_degrees(dt::GW_decorated_graph, e::Edge) = edgeMult(e, dt)
@inline _decorated_valency(dt::GW_decorated_tree, v::Int) = degree(dt.tree, v)
@inline _decorated_valency(dt::GW_decorated_graph, v::Int) = valency(v, dt)

@inline function _bundle_top_weight(E::AbstractGKMVectorBundle, v::Int, t)
  result = one(t[1])
  @inbounds for i in 1:rank(E)
    result *= _bundle_weight_class(E, v, i, t)
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a >= 0 "virtual_zero_section is only implemented for convex vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in 0:(Int(a) * d)
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a < 0 "derivated_functor is only implemented for concave vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in (Int(a) * d + 1):-1
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result = result // _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result *= _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _gw_bundle_functor(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle, mode)
  @req baseof(E) == dt.gkm "The vector bundle and decorated graph must have the same base"
  if dt isa GW_decorated_graph
    @req all(iszero, dt.genus) "Only implemented for genus zero"
  end
  result = _bundle_edge_contribution(one(dt.context.t[1]), dt, E, mode)
  return _bundle_vertex_contribution(result, dt, E, mode)
end

@doc raw"""
    virtual_zero_section(V::AbstractGKMVectorBundle) -> EquivariantClass

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{0,n}(X,\beta)$ of the top Chern class of $\pi_*(\text{ev}^*_{n+1}(V))$ where $\pi\colon \overline{\mathcal{M}}_{0,n+1}(X,\gamma)\rightarrow \overline{\mathcal{M}}_{0,n}(X,\gamma)$ forgets the last map.

!!! note
    This procedure assumes that the moduli space is of stable maps of genus zero and that the vector bundle is convex, i.e., $H^1(\mathbb{P}^1, f^*V) = 0$ for all stable maps $f:\mathbb{P}^1\to X$.
    If these conditions are not met, the computation will stop. If $i\colon Y\hookrightarrow X$ is the zero locus of a generic section of $V$, then the Gromov-Witten invariants of $Y$ satisfy:
    ```math
        \sum_{\gamma:i_*(\gamma)=\beta} i_*([\overline{\mathcal{M}}_{0,n}(Y,\gamma)]^{\mathrm{vir}} = c_{\mathrm{top}}(\pi_*(\text{ev}^*_{n+1}(V)))\cap[\overline{\mathcal{M}}_{0,n}(X,\beta)]^{\mathrm{vir}}
    ```
    See [MR1837116](@cite) and [MR1958379](@cite) for more details.


# Example: Complete intersections in toric varieties
In the following examples, we use **Oscar** to construct a toric variety $X$. At the end, we compute the invariants of $Y$ defined as the zero section of a sum of line bundles of $X$.

Let us start with $X=\mathbb{P}^7$.
```jldoctest P7 
julia> P7 = projective_space(NormalToricVariety, 7)
Normal toric variety
```
The variety $Y$ is the zero section of $V=\mathcal{O}_{\mathbb{P}^7}(2)^{\oplus 4}$. So let us define the Oscar's toric line bundle $L=\mathcal{O}_{\mathbb{P}^7}(2)$.
```jldoctest P7 
julia> L = toric_line_bundle(P7, picard_group(P7)([2]))
Toric line bundle on a normal toric variety
```
Now, we define $l$ to be the **GKM** line bundle defined by $L$. We use [`gkm_line_bundle_of_toric`](@ref).
```jldoctest P7 
julia> l = gkm_line_bundle_of_toric(L)
GKM vector bundle of rank 1 over 8 vertices with weights:
1: (0, 0, 0, 0, 0, 0, 0, 1, -2)
2: (-2, 0, 0, 0, 0, 0, 0, 1, 0)
3: (0, -2, 0, 0, 0, 0, 0, 1, 0)
4: (0, 0, -2, 0, 0, 0, 0, 1, 0)
5: (0, 0, 0, -2, 0, 0, 0, 1, 0)
6: (0, 0, 0, 0, -2, 0, 0, 1, 0)
7: (0, 0, 0, 0, 0, -2, 0, 1, 0)
8: (0, 0, 0, 0, 0, 0, -2, 1, 0)
```
Now we can define the GKM vector bundle $V$.
```jldoctest P7 
julia> V = l + l + l + l;
```
Now we need $\mathbb{P}^7$ defined as a GKM graph.
```jldoctest P7 
julia> X = baseof(V);
```
Finally, we can compute our invariants as usual
```jldoctest P7 
julia> beta = curve_class(X, "1", "2");

julia> P = virtual_zero_section(V);

julia> gromov_witten(X, beta, 0, P; show_bar = false, fast_mode = true)
512

julia> gromov_witten(X, 2*beta, 0, P; show_bar = false, fast_mode = true)
9792
```

Now, let us compute the invariants of $Y$ where $X=\mathbb{P}(\mathcal{O}_{\mathbb{P}^3}\oplus \mathcal{O}_{\mathbb{P}^3}(3))$ and $V=\mathrm{det}(T_X)$.
```jldoctest P3
julia> P3 = projective_space(NormalToricVariety, 3);

julia> O = trivial_line_bundle(P3);

julia> O3 = toric_line_bundle(P3, picard_group(P3)([3]));

julia> Fano = projectivization(O, O3); # Oscar's toric variety X

julia> L = anticanonical_bundle(Fano); # Oscar's anticanonical bundle

julia> V = gkm_line_bundle_of_toric(L) # GKM line bundle
GKM vector bundle of rank 1 over 8 vertices with weights:
1: (0, 0, 0, 0, 1, -2, -1)
2: (0, 0, 0, -2, 1, 0, -7)
3: (-1, 0, 0, 0, 1, -2, 0)
4: (-7, 0, 0, -2, 1, 0, 0)
5: (0, -1, 0, 0, 1, -2, 0)
6: (0, -7, 0, -2, 1, 0, 0)
7: (0, 0, -1, 0, 1, -2, 0)
8: (0, 0, -7, -2, 1, 0, 0)
```

The base of $V$, that is the GKM variety $X$, has effective cone generated by two curve classes.
```jldoctest P3
julia> X = baseof(V);

julia> print_curve_classes(X)
2 -> 1: (0, 1), Chern number: 2
3 -> 1: (1, -3), Chern number: 1
4 -> 2: (1, 0), Chern number: 7
4 -> 3: (0, 1), Chern number: 2
5 -> 1: (1, -3), Chern number: 1
5 -> 3: (1, -3), Chern number: 1
6 -> 2: (1, 0), Chern number: 7
6 -> 4: (1, 0), Chern number: 7
6 -> 5: (0, 1), Chern number: 2
7 -> 1: (1, -3), Chern number: 1
7 -> 3: (1, -3), Chern number: 1
7 -> 5: (1, -3), Chern number: 1
8 -> 2: (1, 0), Chern number: 7
8 -> 4: (1, 0), Chern number: 7
8 -> 6: (1, 0), Chern number: 7
8 -> 7: (0, 1), Chern number: 2

julia> beta_1 = curve_class(X, "1", "2"); # first generator

julia> beta_2 = curve_class(X, "1", "3"); # second generator

julia> P = virtual_zero_section(V);

julia> gromov_witten(X, beta_1, 0, P; show_bar = false, fast_mode = true)
28

julia> gromov_witten(X, beta_2, 0, P; show_bar = false, fast_mode = true)
3
```

# Example: Zero section of homogeneous bundles of Grassmannians
Let us compute the Gromov-Witten invariants of the quintic in $\mathbb{P}^4$. We use [`tautological_and_univ_bd`](@ref).
```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 1, 5);

julia> V = dual(S)^5;

julia> P4 = baseof(V);

julia> beta = curve_class(P4, "1", "2"); # line class

julia> P = virtual_zero_section(V);

julia> gromov_witten(P4, beta, 0, P; show_bar = false, fast_mode = true) # lines in the quintic in P4 
2875

julia> gromov_witten(P4, 2*beta, 0, P; show_bar = false, fast_mode = true)
4876875//8

julia> gromov_witten(P4, 3*beta, 0, P; show_bar = false, fast_mode = true)
8564575000//27
```

Let us compute the Gromov-Witten invariants of the Calabi-Yau threefold given by a zero section of $\mathcal{O}(4)$ in $G(2, 4)$.
```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 2, 4);

julia> plucker = dual(det(S));

julia> V = plucker^4;

julia> G24 = baseof(V);

julia> beta = curve_class(G24, "12", "13");

julia> P = virtual_zero_section(V);

julia> gromov_witten(G24, beta, 0, P; show_bar = false, fast_mode = true)
1280

julia> gromov_witten(G24, 2*beta, 0, P; show_bar = false, fast_mode = true)
92448

julia> gromov_witten(G24, 3*beta, 0, P; show_bar = false, fast_mode = true)
422690816//27
```

Let us compute the Gromov-Witten invariants of the Calabi-Yau threefolds given by a zero section in $G(2, 5)$ of the following bundles:
```math
\begin{aligned}
V1 &= \mathcal{O}(1)\oplus\mathcal{O}(2)\oplus\mathcal{O}(2), \\
V2 &= \mathcal{O}(1)\oplus\mathcal{O}(1)\oplus\mathcal{O}(3), \\
V3 &= \mathcal{S}^{\vee}(1)\oplus\mathcal{O}(2), \\
V4 &= \wedge^2\mathcal{Q}(1).
\end{aligned}
```

```jldoctest
julia> S, Q = tautological_and_univ_bd(GKMGraph, 2, 5);

julia> plucker = dual(det(S));

julia> V1 = plucker + (plucker^2) + (plucker^2);

julia> V2 = plucker + plucker + (plucker^3);

julia> V3 = dual(S) * plucker + (plucker^2);

julia> V4 = wedge_product(Q, 2) * plucker;

julia> P = virtual_zero_section.([V1, V2, V3, V4]);

julia> G25 = baseof(V1);

julia> beta = curve_class(G25, "12", "13");

julia> gromov_witten(G25, beta, 0, P; show_bar = false, fast_mode = true) 
4-element Vector{QQFieldElem}:
 400
 540
 336
 325

julia> gromov_witten(G25, 2*beta, 0, P; show_bar = false, fast_mode = true) 
4-element Vector{QQFieldElem}:
  5590
 25245//2
  3678
 25925//8
```
"""
function virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConvexBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConvexBundleFunctor()))
end

@doc raw"""
    derivated_functor(V::AbstractGKMVectorBundle) -> EquivariantClass

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{0,n}(X,\beta)$ of the top Chern class of $R^{1}\pi_*(\text{ev}^*_{n+1}(V))$ where $\pi\colon \overline{\mathcal{M}}_{0,n+1}(X,\gamma)\rightarrow \overline{\mathcal{M}}_{0,n}(X,\gamma)$ forgets the last map.

!!! note
    This procedure assumes that the moduli space is of stable maps of genus zero and that the vector bundle is concave, i.e., $H^0(\mathbb{P}^1, f^*V) = 0$ for all stable maps $f:\mathbb{P}^1\to X$.
    If these conditions are not met, the computation will stop.
 
# Example: Manin formula
Let us compute the Manin formula, that is:
```math
\int_{\overline{\mathcal{M}}_{0,0}(\mathbb{P}^1, d\beta)} c_{\mathrm{top}}(R^{1}\pi_*(\text{ev}^*_{n+1}(\mathcal{O}_{\mathbb{P}^1}(-1) \oplus \mathcal{O}_{\mathbb{P}^1}(-1)))) = \frac{1}{d^3}.
```
```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 1, 2); # S is the Serre's twisting bundle on P1

julia> V = S + S;

julia> P1 = baseof(V);

julia> beta = curve_class(P1, "1", "2"); # line class

julia> P = derivated_functor(V);

julia> gromov_witten(P1, beta, 0, P; show_bar = false, fast_mode = true)
1

julia> gromov_witten(P1, 2*beta, 0, P; show_bar = false, fast_mode = true)
1//8

julia> gromov_witten(P1, 3*beta, 0, P; show_bar = false, fast_mode = true)
1//27
```
"""
function derivated_functor(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConcaveBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConcaveBundleFunctor()))
end

@doc raw"""
    reduced_virtual_zero_section(V::AbstractGKMVectorBundle) -> EquivariantClass

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{g,n}(X,\beta)$ of the top Chern class of $\pi_*(\text{ev}^*_{n+1}(V))'$, the subbundle of $\pi_*(\text{ev}^*_{n+1}(V))$
that vanishes at the last marked point (cf. [MR1685628; Equation (19)](@cite)).

!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

# Example

Let us use this function to compute the equivariant twisted quantum product $[p_1]\ast_{\mathcal{O}(2)} [p_2]$ on $\mathbb{P}^3$ with respect to the convex vector bundle $\mathcal{O}(2)$, as defined in [MR1685628; Seciton 2.1](@cite).
We compute the coefficients up to degree 3.

```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3);

julia> L = toric_line_bundle(P3, picard_group(P3)([2]));

julia> l = gkm_line_bundle_of_toric(L);

julia> X = baseof(l);

julia> p1 = point_class(X, 1)
GKM class with restrictions: 
[t1*t2*t3 - t1*t2*t5 - t1*t3*t5 + t1*t5^2 - t2*t3*t5 + t2*t5^2 + t3*t5^2 - t5^3, 0, 0, 0]

julia> p2 = point_class(X, 2)
GKM class with restrictions: 
[0, -t1^3 + t1^2*t2 + t1^2*t3 + t1^2*t5 - t1*t2*t3 - t1*t2*t5 - t1*t3*t5 + t2*t3*t5, 0, 0]

julia> beta = curve_class(X, "1", "2");

julia> P = reduced_virtual_zero_section(l);

julia> quantum_product(X, beta, p1, p2; twist_class=P)
GKM class with restrictions: 
[2*t1^2*t2*t3 - 2*t1^2*t2*t5 - 2*t1^2*t3*t5 + 2*t1^2*t5^2 - 3*t1*t2*t3*t4 + 2*t1*t2*t3*t5 + 3*t1*t2*t4*t5 - 2*t1*t2*t5^2 + 3*t1*t3*t4*t5 - 2*t1*t3*t5^2 - 3*t1*t4*t5^2 + 2*t1*t5^3 + t2*t3*t4^2 - t2*t3*t4*t5 - t2*t4^2*t5 + t2*t4*t5^2 - t3*t4^2*t5 + t3*t4*t5^2 + t4^2*t5^2 - t4*t5^3, -t1^3*t4 + 2*t1^3*t5 + t1^2*t2*t4 - 2*t1^2*t2*t5 + t1^2*t3*t4 - 2*t1^2*t3*t5 + t1^2*t4^2 - 3*t1^2*t4*t5 + 2*t1^2*t5^2 - t1*t2*t3*t4 + 2*t1*t2*t3*t5 - t1*t2*t4^2 + 3*t1*t2*t4*t5 - 2*t1*t2*t5^2 - t1*t3*t4^2 + 3*t1*t3*t4*t5 - 2*t1*t3*t5^2 + t2*t3*t4^2 - 3*t2*t3*t4*t5 + 2*t2*t3*t5^2, 0, 0]

julia> quantum_product(X, 2*beta, p1, p2; twist_class=P)
GKM class with restrictions: 
[-2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5]

julia> quantum_product(X, 3*beta, p1, p2; twist_class=P)
GKM class with restrictions: 
[0, 0, 0, 0]
```
Note that one of the equivariant parameters in this examples comes from the fibrewise scaling of $\mathcal{O}(2)$.
```jldoctest
julia> V = vector_bundle_O(2, [1]);

julia> P2 = baseof(V);

julia> beta = curve_class(P2, Edge(1, 2));

julia> P = ev(1, point_class(P2, 1)) * reduced_virtual_zero_section(V);

julia> gromov_witten(P2, beta, 1, P; fast_mode = true, show_bar = false)
1
```
"""
function reduced_virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_reduced_virtual_zero_section(dt, $E))
  return EquivariantClass(rule, dt -> _reduced_virtual_zero_section(dt, E))
end

function _reduced_virtual_zero_section(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle)
  @req !isempty(dt.marks) "Need at least one marked point to reduce the virtual zero section"
  result = _gw_bundle_functor(dt, E, _ConvexBundleFunctor())
  marked_vertex = imageOf(last(dt.marks), dt)
  return result // _bundle_top_weight(E, marked_vertex, dt.context.t)
end
