# Positive genus version of /src/GW/Euler.jl.
function Euler_inv_pos_gen(dg::GW_decorated_graph, t::Vector{T}, edge_weight_dict::Dict{Edge, T}, point_weight_dict::Vector{Union{Nothing, T}}, H::Dict{HodgeKey, QQFieldElem}) where T<:RingElem

  res = one(t[1])

  for v in 1:n_vertices(dg.g)

    valv_plus_g = degree(dg.g, v) + dg.genus[v]
    e = euler_class(imageOf(v, dg), dg.gkm.equivariantCohomology, t, edge_weight_dict, point_weight_dict)
    #println("e = $e, val = $valv")
    if valv_plus_g >= 1
      res = res * e^(valv_plus_g - 1)
    else
      res = res // e
    end

    ### EXPERIMENTAL PART - not justified by the Liu--Sheshmani formula:
    # res = res // prod(edgeMult(Edge(v,n), dg) for n in all_neighbors(dg.g, v))^(2*dg.genus[v])
    ### END OF EXPERIMENTAL PART.

    imV = imageOf(v, dg)
    u = [edgeMult(Edge(v, n), dg) // weight_class(Edge(imV, imageOf(n, dg)), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.g, v)]
    w = [1 // weight_class(Edge(imV, n), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.gkm.g, imV)]
    
    nMarks = count(i -> i==v, dg.marks)
    vpEval = evaluate_vertex_polynomial(u, w, nMarks, dg.genus[v], H)
    res = res * vpEval

    if dg.genus[v] > 0
     println("vpEval = $(factor(numerator(vpEval))) // $(factor(denominator(vpEval)))")
     println("u = $u")
     println("w = $w")
     println("nMarks = $nMarks")
    end
  end

  # if n_edges(dg.g) == 1 && sum(dg.genus) == 1
  #   res = res // 4
  # end

  return res
end

function Euler_inv_pos_gen_VB(V::GKM_vector_bundle, dg::GW_decorated_graph, H::Dict{HodgeKey, QQFieldElem})

  R = V.gkm.equivariantCohomology.coeffRing

  res = one(R)

  for v in 1:n_vertices(dg.g)

    valv_plus_g = degree(dg.g, v) + dg.genus[v]
    e = euler_class(imageOf(v, dg), dg.gkm.equivariantCohomology)
    e = e * _fiber_normal_weight(imageOf(v, dg), V)
    #println("e = $e, val = $valv")
    if valv_plus_g >= 1
      res = res * e^(valv_plus_g - 1)
    else
      res = res // e
    end

    ### EXPERIMENTAL PART - not justified by the Liu--Sheshmani formula:
    # res = res // prod(edgeMult(Edge(v,n), dg) for n in all_neighbors(dg.g, v))^(2*dg.genus[v])
    ### END OF EXPERIMENTAL PART.

    imV = imageOf(v, dg)
    u = [edgeMult(Edge(v, n), dg) // weight_class(Edge(imV, imageOf(n, dg)), dg.gkm) for n in all_neighbors(dg.g, v)]
    w_from_base = [1 // weight_class(Edge(imV, n), dg.gkm) for n in all_neighbors(dg.gkm.g, imV)]
    w_from_fibre = [1 // _fiber_summand_weight(imV, i, V) for i in 1:rank(V)]
    w = vcat(w_from_base, w_from_fibre)
    
    nMarks = count(i -> i==v, dg.marks)
    vpEval = evaluate_vertex_polynomial(u, w, nMarks, dg.genus[v], H)
    res = res * vpEval

    # if dg.genus[v] > 0
    #  println("vpEval = $(factor(numerator(vpEval))) // $(factor(denominator(vpEval)))")
    #  println("u = $u")
    #  println("w = $w")
    #  println("nMarks = $nMarks")
    # end
  end

  # if n_edges(dg.g) == 1 && sum(dg.genus) == 1
  #   res = res // 4
  # end

  return res
end