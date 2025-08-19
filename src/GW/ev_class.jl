
@doc raw"""
    ev(j::Int64, cc) -> EquivariantClass

Return the equivariant cohomology class on $\overline{\mathcal{M}_{0,n}}(X,\beta)$ that is given by pulling back
the cohomology class `cc` on $X$ along the evaluation map at the `j`-th point: $\text{ev}_j\colon \overline{\mathcal{M}_{0,n}}(X,\beta)\rightarrow X$.

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