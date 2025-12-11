export virtual_zero_section
export reduced_virtual_zero_section

@doc raw"""
    virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}}_{g,n}(X,\beta)$ of the top Chern class of $\pi_*(\text{ev}^*_{n+1}(V))$ where $\pi\colon \overline{\mathcal{M}}_{0,n+1}(X,\gamma)\rightarrow \overline{\mathcal{M}}_{0,n}(X,\gamma)$ forgets the last map.

!!! note
    This procedure assumes that the moduli space is of stable maps of genus zero and that the vector bundle is convex, i.e., $H^1(\mathbb{P}^1, f^*V) = 0$ for all stable maps $f:\mathbb{P}^1\to X$.
    If these conditions are not met, the computation will stop. If $i\colon Y\hookrightarrow X$ is the zero locus of a generic section of $V$, then the Gromov-Witten invariants of $Y$ satisfy:
    ```math
        \sum_{\gamma:i_*(\gamma)=\beta} i_*([\overline{\mathcal{M}}_{0,n}(Y,\gamma)]^{\mathrm{vir}} = c_{\mathrm{top}}(\pi_*(\text{ev}^*_{n+1}(V)))\cap[\overline{\mathcal{M}}_{0,n}(X,\beta)]^{\mathrm{vir}}
    ```
    See [MR1837116](@cite) and [MR1958379](@cite) for more details.


# Example
Let us compute the Gromov-Witten invariants of the quintic in P4. We use [`tautological_and_univ_bd`](@ref).
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

julia> beta = curve_class(G24, "12", "13"); # line class

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

  rule = :(_virtual_zero_section(dt, $V))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _virtual_zero_section(dt::Union{GW_decorated_tree, GW_decorated_graph}, V::GKM_vector_bundle)

  @req V.gkm == dt.gkm "The vector bundle V must be over the same GKM graph as dt.gkm."

  if isa(dt, GW_decorated_graph)
    if any(i-> dt.genus[i] > 0, eachindex(dt.genus))
      error("virtual_zero_section only implemented for genus zero.")
    end
  end

  C = get_connection(V)
  @req !isnothing(C) "Vector bundle needs a connection" #TODO: this could be any compatible connection.
  _calculate_connection_a(V)
  _calculate_weight_classes(V)

  R = dt.gkm.equivariantCohomology.coeffRing
  rV = rank(V)
  ans = one(R)

  for e in edges(dt.tree)
    
    # color endpoints
    v1 = imageOf(src(e), dt)
    v2 = imageOf(dst(e), dt)

    # compute vector [a1, a2, ..., ar]
    a_vector = [_fiber_connection_a(imageOf(e, dt), i, V) for i in 1:rV]

    @req all(i -> a_vector[i] >= 0, eachindex(a_vector)) "virtual_zero_section only implemented for convex vector bundles, got $a_vector."
    # The error message above is not equivalent to the condition that is checked!
    # The error message asks for a1+...+ar > 0, while the condition checks a1>0 && ... && ar>0.

    de = dt.edgeMult[e]
    lambda = (sum(i -> V.gkm.w[Edge(v1,v2)][i]*gens(R)[i], eachindex(gens(R))))//de; # weight of the tangent bundle, divided by de

    for i in 1:rV
      a = a_vector[i]

      lambda_1 = _fiber_summand_weight(v1, i, V)
      # lambda_2 = _fiber_summand_weight(v2, C[(e ,i)], V) 

      # original formulation, work for line bundles
      # for k in 0:(Int(a)*de)
      #   ans *= (k*lambda_1 + (a*de-k) * lambda_2) // (a*de) ; push!(vec_we, (k*lambda_1 + (a*de-k) * lambda_2))
      # end

      # new formulation
      for k in 0:(Int(a)*de)
        ans *= lambda_1 - k*lambda
      end
    end
    
  end

  for v in 1:n_vertices(dt.tree)
    val = length(all_neighbors(dt.tree, v))
    val == 1 && continue
  
    # with the following two lines it works for line bundles
    # lambda_v = gens(R)[imageOf(v, dt)]
    # ans *= (1//(5*lambda_v))^(-(1 - val))
    #############################

    lambda_v = _fiber_normal_weight(imageOf(v, dt), V)#; println(lambda_v)
    ans *= (1//(lambda_v))^(-(1 - val))  
  end
    
  return ans
end


@doc raw"""
    reduced_virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}_{g,n}}(X,\beta)$ of the top Chern class of the subbundle of $\pi_*(\text{ev}^*_{n+1}(V))$
that vanishes at the last marked point (cf. [MR1685628; Equation (19)](@cite)).

!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.
"""
function reduced_virtual_zero_section(V::GKM_vector_bundle)::EquivariantClass

  rule = :(_reduced_virtual_zero_section(dt, $V))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

# Like _victual_zero_section, but divide out by the top chern class of the vector bundle at the last marked point.
function _reduced_virtual_zero_section(dt::Union{GW_decorated_tree, GW_decorated_graph}, V::GKM_vector_bundle)

  ans = _virtual_zero_section(dt, V)

  @req length(dt.marks) >= 1 "Need at least one marked point to reduce the virtual zero section."  

  ans //= _fiber_normal_weight(imageOf(dt.marks[length(dt.marks)], dt), V)

  return ans
end
