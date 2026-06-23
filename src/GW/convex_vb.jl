export virtual_zero_section
export reduced_virtual_zero_section
export derivated_functor

##############################
# Ultra-optimized GW functors (Julia)
# Allocation-free, dispatch-based, hoisted invariants
##############################

# ==========================================================
# Mode tags (compile-time dispatch)
# ==========================================================
struct Convex end
struct Concave end

# ==========================================================
# Internal helpers (fully inlined)
# ==========================================================

# TODO: write down the formula that this computes and compare it to Overleaf.
@inline function _tangent_weight(dt, V, v1, v2, de)
    R     = dt.gkm.equivariantCohomology.coeffRing
    gensR = gens(R)
    w     = _w(V.gkm, Edge(v1, v2))
    λ = zero(R)
    @inbounds @simd for i in eachindex(gensR)
        λ += w[i] * gensR[i]
    end
    return λ // de
end

# ==========================================================
# Edge contributions (hot kernels)
# ==========================================================

function _edge_contribution!(ans, dt, V, ::Convex)
    rV = rank(V)

    for e in edges(dt.tree)
        v1 = imageOf(src(e), dt)
        v2 = imageOf(dst(e), dt)
        de = dt.edgeMult[e]

        λ = _tangent_weight(dt, V, v1, v2, de)

        @inbounds for i in 1:rV
            a = _fiber_connection_a(imageOf(e, dt), i, V)
            @req a ≥ 0 "virtual_zero_section only implemented for convex vector bundles."

            λ1 = _fiber_summand_weight(v1, i, V)
            N  = Int(a) * de

            @inbounds @simd for k in 0:N
                ans *= λ1 - k * λ
            end
        end
    end
    return ans
end


function _edge_contribution!(ans, dt, V, ::Concave)
    rV = rank(V)

    for e in edges(dt.tree)
        v1 = imageOf(src(e), dt)
        v2 = imageOf(dst(e), dt)
        de = dt.edgeMult[e]

        λ = _tangent_weight(dt, V, v1, v2, de)

        @inbounds for i in 1:rV
            a = _fiber_connection_a(imageOf(e, dt), i, V)
            @req a < 0 "derivated_functor only implemented for concave vector bundles."

            λ1 = _fiber_summand_weight(v1, i, V)
            N  = Int(a) * de

            @inbounds @simd for k in (N + 1):-1
                ans *= λ1 - k * λ
            end
        end
    end
    return ans
end

# ==========================================================
# Vertex contributions
# ==========================================================

function _vertex_contribution(ans, dt, V, ::Convex)
    for v in 1:n_vertices(dt.tree)
        val = length(all_neighbors(dt.tree, v))
        val == 1 && continue
        λv = _fiber_normal_weight(imageOf(v, dt), V)
        ans *= (1 // λv)^(val - 1)
    end
    return ans
end


function _vertex_contribution(ans, dt, V, ::Concave)
    for v in 1:n_vertices(dt.tree)
        val = length(all_neighbors(dt.tree, v))
        val == 1 && continue
        λv = _fiber_normal_weight(imageOf(v, dt), V)
        ans *= λv^(val - 1)
    end
    return ans
end

# ==========================================================
# Core mule function (single specialization point)
# ==========================================================

function _gw_functor_core(
    dt::Union{GW_decorated_tree, GW_decorated_graph},
    V::GKM_vector_bundle,
    mode
)
    @req V.gkm == dt.gkm "The vector bundle V must be over the same GKM graph as dt.gkm."

    if isa(dt, GW_decorated_graph)
        any(i -> dt.genus[i] > 0, eachindex(dt.genus)) &&
            error("Only implemented for genus zero.")
    end


    R   = dt.gkm.equivariantCohomology.coeffRing
    ans = one(R)

    ans = _edge_contribution!(ans, dt, V, mode)
    ans = _vertex_contribution(ans, dt, V, mode)

    return ans
end

# ==========================================================
# Public API (zero-cost wrappers)
# ==========================================================

# @inline _virtual_zero_section(dt, V) =
#     _gw_functor_core(dt, V, Convex())

# @inline _derivated_functor(dt, V) =
#     _gw_functor_core(dt, V, Concave())

function _connection_ready(V::GKM_vector_bundle)
  
  C = get_connection(V)
  @req !isnothing(C) "Vector bundle needs a connection"

  _calculate_connection_a(V)
  _calculate_weight_classes(V)
  return
end

@doc raw"""
    virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

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
GKM vector bundle of rank 1 over GKM graph with 8 nodes and valency 7 with weights:
1: (0, 0, 0, 0, 0, 0, 0, -1, 2)
2: (2, 0, 0, 0, 0, 0, 0, -1, 0)
3: (0, 2, 0, 0, 0, 0, 0, -1, 0)
4: (0, 0, 2, 0, 0, 0, 0, -1, 0)
5: (0, 0, 0, 2, 0, 0, 0, -1, 0)
6: (0, 0, 0, 0, 2, 0, 0, -1, 0)
7: (0, 0, 0, 0, 0, 2, 0, -1, 0)
8: (0, 0, 0, 0, 0, 0, 2, -1, 0)
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
GKM vector bundle of rank 1 over GKM graph with 8 nodes and valency 4 with weights:
1: (0, 0, 0, 0, -1, 2, 1)
2: (0, 0, 0, 2, -1, 0, 7)
3: (1, 0, 0, 0, -1, 2, 0)
4: (7, 0, 0, 2, -1, 0, 0)
5: (0, 1, 0, 0, -1, 2, 0)
6: (0, 7, 0, 2, -1, 0, 0)
7: (0, 0, 1, 0, -1, 2, 0)
8: (0, 0, 7, 2, -1, 0, 0)
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
julia> S, _ = tautological_and_univ_bd(GKM_graph, 1, 5);

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
julia> S, _ = tautological_and_univ_bd(GKM_graph, 2, 4);

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
julia> S, Q = tautological_and_univ_bd(GKM_graph, 2, 5);

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
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function virtual_zero_section(V::GKM_vector_bundle)::EquivariantClass

  _connection_ready(V)

  rule = :(_gw_functor_core(dt, $V, Convex()))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

@doc raw"""
    derivated_functor(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{0,n}(X,\beta)$ of the top Chern class of $R^{1}\pi_*(\text{ev}^*_{n+1}(V))$ where $\pi\colon \overline{\mathcal{M}}_{0,n+1}(X,\gamma)\rightarrow \overline{\mathcal{M}}_{0,n}(X,\gamma)$ forgets the last map.

!!! note
    This procedure assumes that the moduli space is of stable maps of genus zero and that the vector bundle is concave, i.e., $H^0(\mathbb{P}^1, f^*V) = 0$ for all stable maps $f:\mathbb{P}^1\to X$.
    If these conditions are not met, the computation will stop.
 
# Example: Manin formula
Let us compute the Manin formula, that is:
```math
\int_{\overline{\mathcal{M}}_{0,0}(\mathbb{P}^1, d\beta)} R^{1}\pi_*(\text{ev}^*_{n+1}(\mathcal{O}_{\mathbb{P}^1}(-1) \oplus \mathcal{O}_{\mathbb{P}^1}(-1))) = \frac{1}{d^3}.
```
```jldoctest
julia> S, _ = tautological_and_univ_bd(GKM_graph, 1, 2); # S is the Serre's twisting bundle on P1

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
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.
"""
function derivated_functor(V::GKM_vector_bundle)::EquivariantClass
  
  _connection_ready(V)
  
  rule = :(_gw_functor_core(dt, $V, Concave()))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

@doc raw"""
    reduced_virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A _convex_ vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{g,n}(X,\beta)$ of the top Chern class of the subbundle of $\pi_*(\text{ev}^*_{n+1}(V))$
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
(-t1*t2*t3 + t1*t2*t5 + t1*t3*t5 - t1*t5^2 + t2*t3*t5 - t2*t5^2 - t3*t5^2 + t5^3)*e[1]

julia> p2 = point_class(X, 2)
(t1^3 - t1^2*t2 - t1^2*t3 - t1^2*t5 + t1*t2*t3 + t1*t2*t5 + t1*t3*t5 - t2*t3*t5)*e[2]

julia> beta = curve_class(X, "1", "2");

julia> P = reduced_virtual_zero_section(l)

julia> quantum_product(X, beta, p1, p2; useStructureConstants=false, twist_class=P)
(2*t1^2*t2*t3 - 2*t1^2*t2*t5 - 2*t1^2*t3*t5 + 2*t1^2*t5^2 - 3*t1*t2*t3*t4 + 2*t1*t2*t3*t5 + 3*t1*t2*t4*t5 - 2*t1*t2*t5^2 + 3*t1*t3*t4*t5 - 2*t1*t3*t5^2 - 3*t1*t4*t5^2 + 2*t1*t5^3 + t2*t3*t4^2 - t2*t3*t4*t5 - t2*t4^2*t5 + t2*t4*t5^2 - t3*t4^2*t5 + t3*t4*t5^2 + t4^2*t5^2 - t4*t5^3, -t1^3*t4 + 2*t1^3*t5 + t1^2*t2*t4 - 2*t1^2*t2*t5 + t1^2*t3*t4 - 2*t1^2*t3*t5 + t1^2*t4^2 - 3*t1^2*t4*t5 + 2*t1^2*t5^2 - t1*t2*t3*t4 + 2*t1*t2*t3*t5 - t1*t2*t4^2 + 3*t1*t2*t4*t5 - 2*t1*t2*t5^2 - t1*t3*t4^2 + 3*t1*t3*t4*t5 - 2*t1*t3*t5^2 + t2*t3*t4^2 - 3*t2*t3*t4*t5 + 2*t2*t3*t5^2, 0, 0)

julia> quantum_product(X, 2*beta, p1, p2; useStructureConstants=false, twist_class=P)
(-2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5, -2*t1*t4 + 4*t1*t5 + t4^2 - 2*t4*t5)

julia> quantum_product(X, 3*beta, p1, p2; useStructureConstants=false, twist_class=P)
(0, 0, 0, 0)
```
Note that one of the equivariant parameters in this examples comes from the fibrewise scaling of $\mathcal{O}(2)$.
"""
function reduced_virtual_zero_section(V::GKM_vector_bundle)::EquivariantClass

    _connection_ready(V)

  rule = :(_reduced_virtual_zero_section(dt, $V))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

# Like _victual_zero_section, but divide out by the top chern class of the vector bundle at the last marked point.
function _reduced_virtual_zero_section(dt::Union{GW_decorated_tree, GW_decorated_graph}, V::GKM_vector_bundle)

  ans = _gw_functor_core(dt, V, Convex())

  @req length(dt.marks) >= 1 "Need at least one marked point to reduce the virtual zero section."  

  ans //= _fiber_normal_weight(imageOf(dt.marks[length(dt.marks)], dt), V)

  return ans
end