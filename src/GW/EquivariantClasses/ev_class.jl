@doc raw"""
    ev(j::Int64, cc) -> EquivariantClass

Return the equivariant cohomology class on $\overline{\mathcal{M}_{0,n}}(X,\beta)$ that is given by pulling back
the cohomology class `cc` on $X$ along the evaluation map at the `j`-th point: $\text{ev}_j\colon \overline{\mathcal{M}_{0,n}}(X,\beta)\rightarrow X$.

For example, if `cc` is a point class, then `ev(j, cc)` is the pullback of the point class at the `j`-th marked point. The following invariant:

```math
\int_{\overline{M}_{0,2}(\mathbb{P}^3, 1)}\mathrm{ev}_{1}^{*}([\mathrm{pt}])\cdot\mathrm{ev}_{2}^{*}([\mathrm{pt}]) = 1,
```

can be computed as follows:
```jldoctest
julia> G = projective_space(GKM_graph, 3);

julia> P = ev(1, point_class(G, 1)) * ev(2, point_class(G, 2));

julia> beta = curve_class(G, "1", "2"); # beta is the class of a line in $\mathbb{P}^3$

julia> gromov_witten(G, beta, 2, P; show_bar=false)
1
```

The following is a more complex example, where we compute the Gromov-Witten invariant of the Grassmannian $G(2,4)$.
# Example
```jldoctest ev
julia> G24 = grassmannian(GKM_graph, 2, 4);

julia> e1 = ev(1, point_class(G24, 1));

julia> e2 = ev(2, first_chern_class(G24));

julia> e3 = ev(3, poincare_dual(gkm_subgraph_from_vertices(G24, [1, 2])));

julia> beta = curve_class(G24, Edge(1, 2));

julia> gromov_witten(G24, beta, 3, e1 * e2 * e3; show_bar=false)
4
```
"""
function ev(j::Int64, cc)::EquivariantClass

  rule = :(_ev(dt, $j, $cc))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _ev(dt::GW_decorated_tree, j::Int64, cc)

  v = imageOf(dt.marks[j], dt)
  
  return cc[v]
end