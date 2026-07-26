function Euler_inv_pos_gen(dg::GW_decorated_graph, t::Vector{T}, psi_exp::Vector{Int}, H::Dict{HodgeKey,QQFieldElem}) where T
  result = one(t[1])
  target_valency = valency(dg.gkm)
  for v in vertices(dg.g)
    genus = dg.genus[v]
    valence = valency(v, dg)
    euler = _euler_class(dg.gkm, imageOf(v, dg), t)
    exponent = valence + genus - 1
    result = exponent >= 0 ? result * euler^exponent : result // euler^(-exponent)

    image_vertex = imageOf(v, dg)
    u = reduce(vcat, (edgeMult(Edge(v, n), dg) .// _weight_class(dg.gkm, Edge(image_vertex, imageOf(n, dg)), t) for n in all_neighbors(dg.g, v)); init=Vector{typeof(one(t[1]) // one(t[1]))}())
    w = [one(t[1]) // _flag_weight_class(dg.gkm, image_vertex, i, t) for i in 1:target_valency]
    psi_at_vertex = Int[psi_exp[i] for i in eachindex(dg.marks) if dg.marks[i] == v]
    result *= evaluate_vertex_polynomial_with_psis(u, w, psi_at_vertex, genus, H)
  end
  return result
end
