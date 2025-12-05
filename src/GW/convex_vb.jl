export virtual_zero_section
export reduced_virtual_zero_section

@doc raw"""
    virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}_{g,n}}(X,\beta)$ of the top Chern class of $\pi_*(\text{ev}^*_{n+1}(V))$ . That is:

```math
c_{\text{top}}(\pi_*(\text{ev}^*_{n+1}(V))) = 
\prod_{e=\{p,q\}\in E(\Gamma)}\left(\prod_{k=0}^{V_e}\frac{k\,c_1(V)|_{\overrightarrow{f}(p)}+(V_e-k)\,c_1(V)|_{\overrightarrow{f}(q)}}{V_e}\right)\prod_{v\in V(\Gamma)}\left(c_1(V)|_{\overrightarrow{f}(v)}\right)^{1-\mathrm{val}(v)},
```
where $\Gamma$ is the decorated graph corresponding to a fixed locus in $\overline{\mathcal{M}_{g,n}}(X,\beta)$, $V_e$ is the edge multiplicity of the edge $e$ in $\Gamma$ times the degree of the intersection $c_1(V)\cup C_e$, and $\mathrm{val}(v)$ is the valency of the vertex $v$ in $\Gamma$.

!!! warning
    This procedure assumes that the moduli space is of stable maps of genus zero and that the vector bundle is convex, i.e., $H^1(\mathbb{P}^1, f^*V) = 0$ for all stable maps $f:\mathbb{P}^1\to X$.
    If these conditions are not met, the result will be incorrect. If $i\colon Y\hookrightarrow X$ is the zero locus of a generic section of $V$, then the Gromov-Witten invariants of $Y$ satisfy:
    ```math
        \sum_{\gamma:i_*(\gamma)=\beta} i_*([\overline{M}_{0,n}(Y,\gamma)]^{\mathrm{vir}} = c_{\mathrm{top}}(\mathbf{V})\cap[\overline{M}_{0,n}(X,\beta)]^{\mathrm{vir}}
    ```
    See [MR1837116](@cite) and [MR1958379](@cite) for more details.


# Example
Let us compute the Gromov-Witten invariants of the quintic in $\mathbb{P}^4$.
```jldoctest
julia> P4 = projective_space(GKM_graph, 4);

julia> beta = curve_class(G24, "1", "2"); # line class

julia> V = 

julia> gromov_witten(G24, beta, 3, e1 * e2 * e3; show_bar=false)
4
```

Note that if $X$ is a toric variety, one may use Oscar's `ToricVectorBundle` to define the vector bundle. Let us compute the Gromov-Witten invariant of the cubic in $\mathbb{P}^3$ using a toric vector bundle.
```jldoctest ev
julia> P3 = projective_space(NormalToricVariety, 3); # P3 given as Oscar's toric variety

julia> l = toric_line_bundle(P3, [ZZRingElem(3)]); # cubic line bundle

julia> X, V = convert_toric_vector_bundle_to_GKM_vector_bundle(l);

julia> beta = curve_class(X, "1", "2"); # line class

julia> gromov_witten(X, beta, 3, e1 * e2 * e3; show_bar=false)
4
```
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

    @req all(i -> a_vector[i] >= 0, eachindex(a_vector)) "virtual_zero_section only implemented for convex vector bundles."
    @req all(i -> a_vector[i] > 0, eachindex(a_vector)) "The first Chern class of the vector bundle must pair positively with curve classes." #TODO: fix formulation.
    # The error message above is not equivalent to the condition that is checked!
    # The error message asks for a1+...+ar > 0, while the condition checks a1>0 && ... && ar>0.
    
    de = dt.edgeMult[e]

    for i in 1:rV
      a = a_vector[i]
      # for k in 0:(Int(a)*de)
      #   ans *= (k * lambda_1 + (a*de - k) * lambda_2) // de
      # end

      lambda_1 = _fiber_summand_weight(v1, i, V)
      lambda_2 = _fiber_summand_weight(v2, C[(e ,i)], V)
      for k in 0:(Int(a)*de)
        ans *= (k*lambda_1 + (a*de-k) * lambda_2) // (a*de)
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

    lambda_v = _fiber_normal_weight(imageOf(v, dt), V); #println(lambda_v)
    ans *= (1//(lambda_v))^(-(1 - val))

    # lambda_v = prod(i -> _fiber_summand_weight(v, i, V), 1:rank(V))
    # ans *= (1//(lambda_v))^(-(1 - val))
  
  end
    
  return ans
end


@doc raw"""
    virtual_zero_section(V::GKM_vector_bundle) -> EquivariantClass

# Arguments
 - `V::GKM_vector_bundle`: A vector bundle over a GKM graph $X$.

Return the equivariant cohomology class on $\overline{\mathcal{M}_{g,n}}(X,\beta)$ of the top Chern class of the subbundle of $\pi_*(\text{ev}^*_{n+1}(V))$
that vanishes at the last marked point (cf. [MR1685628; Equation (19)](@cite)).

# Example
TODO: write down example.
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
