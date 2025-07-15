@doc raw"""
  class_one() -> EquivariantClass

It returns the cohomology class $1$ on $\overline{\mathcal{M}_{0,n}}(X,\beta)$. This is useful when $\overline{\mathcal{M}_{0,n}}(X,\beta)$ is zero dimensional.

# Example
```jldoctest class_one
julia> F = hirzebruch_surface(NormalToricVariety, 1);

julia> G = gkm_graph_of_toric(F); 

julia> beta = curve_class(G, "1", "2"); # class of the zero section of the Hirzebruch surface

julia> gromov_witten(G, beta, 0, class_one(), show_bar = false)
1
```
"""
function class_one()::EquivariantClass

  rule = :(_class_one(dt))

  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _class_one(dt::GW_decorated_tree)

  return 1
end