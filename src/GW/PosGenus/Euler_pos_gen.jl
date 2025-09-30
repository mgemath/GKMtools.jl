# Positive genus version of /src/GW/Euler.jl.
function Euler_inv_pos_gen(dg::GW_decorated_graph, t::Vector{T}, edge_weight_dict::Dict{Edge, T}, point_weight_dict::Vector{Union{Nothing, T}}, H::Dict{HodgeKey, QQFieldElem}) where T<:RingElem

  res = one(t[1])

  for v in 1:n_vertices(dg.g)

    valv = degree(dg.g, v)
    e = euler_class(imageOf(v, dg), dg.gkm.equivariantCohomology, t, edge_weight_dict, point_weight_dict)
    #println("e = $e, val = $valv")
    if valv >= 1
      res = res * e^(valv - 1)
    else
      res = res // e
    end

    imV = imageOf(v, dg)
    u = [edgeMult(Edge(v, n), dg) // weight_class(Edge(imV, imageOf(n, dg)), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.g, v)]
    w = [1 // weight_class(Edge(imV, n), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.gkm.g, imV)]
    nMarks = count(i -> i==v, dg.marks)
    res = res * evaluate_vertex_polynomial(u, w, nMarks, dg.genus[v], H)

  end

  return res
end