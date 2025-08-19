@doc raw"""
    class_one() -> EquivariantClass

Return the cohomology class $1$ on $\overline{\mathcal{M}_{0,n}}(X,\beta)$.

# Example
```jldoctest class_one
julia> P2 = projective_space(GKM_graph, 2);

julia> beta = curve_class(P2, Edge(1, 2));

julia> gromov_witten(P2, beta, 2, ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 2)); show_bar=false)
1

julia> gromov_witten(P2, beta, 2, ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 2)) * class_one(); show_bar=false)
1

julia> gromov_witten(P2, beta, 0, class_one(); show_bar=false)
0
```
"""
function class_one()::EquivariantClass

    rule = :(_class_one(dt))

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _class_one(dt::GW_decorated_tree)

    return 1
end