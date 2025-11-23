# Positive genus version of /src/GW/Euler.jl.
function Euler_inv_pos_gen(dg::GW_decorated_graph, t::Vector{T}, edge_weight_dict::Dict{Edge, T}, point_weight_dict::Vector{Union{Nothing, T}}, VPs::Matrix{QQMPolyRingElem}) where T<:RingElem

  res = one(t[1])

  for v in 1:n_vertices(dg.g)

    valv_plus_g = degree(dg.g, v) + dg.genus[v]
    # valv_plus_g_minus_one = valency(v, dg) + dg.genus[v] - 1
    e = euler_class(imageOf(v, dg), dg.gkm.equivariantCohomology, t, edge_weight_dict, point_weight_dict)

    # res = res * (e^valv_plus_g_minus_one) # This probably throws an error if valv_plus_g_minus_one = -1. That's why I used the six lines below originally.

    #println("e = $e, val = $valv")
    if valv_plus_g >= 1
      res = res * e^(valv_plus_g - 1)
    else
      res = res // e
    end


    imV = imageOf(v, dg)
    u = vcat((edgeMult(Edge(v, n), dg) .// weight_class(Edge(imV, imageOf(n, dg)), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.g, v))...,)
    w = [1 // weight_class(Edge(imV, n), dg.gkm, t, edge_weight_dict) for n in all_neighbors(dg.gkm.g, imV)]
    
    nMarks = count(i -> i==v, dg.marks)
    vpEval = evaluate_vertex_polynomial(u, w, nMarks, dg.genus[v], VPs)
    res = res * vpEval

    # if dg.genus[v] > 0
    #   println("vpEval = $(factor(numerator(vpEval))) // $(factor(denominator(vpEval)))")
    #   println("u = $u")
    #   println("w = $w")
    #   println("nMarks = $nMarks")
    # end
  end

  # if n_edges(dg.g) == 1 && sum(dg.genus) == 1
  #   res = res // 4
  # end

  return res
end

# This returns the extra factor for Euler_inv_pos_gen in the fiber direction.
function _Euler_inv_pos_gen_VB(dg::GW_decorated_graph, V::GKM_vector_bundle, VPs::Matrix{QQMPolyRingElem})::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  R = V.gkm.equivariantCohomology.coeffRing

  res = one(R)

  for v in 1:n_vertices(dg.g)

    valv_plus_g = valency(v, dg) + dg.genus[v]
    e = euler_class(imageOf(v, dg), dg.gkm.equivariantCohomology)
    e = e * _fiber_normal_weight(imageOf(v, dg), V)
    #println("e = $e, val = $valv")
    if valv_plus_g >= 1
      res = res * e^(valv_plus_g - 1)
    else
      res = res // e
    end

    imV = imageOf(v, dg)
    u = vcat((edgeMult(Edge(v, n), dg) .// weight_class(Edge(imV, imageOf(n, dg)), dg.gkm) for n in all_neighbors(dg.g, v))...,)
    w_from_base = [1 // weight_class(Edge(imV, n), dg.gkm) for n in all_neighbors(dg.gkm.g, imV)]
    w_from_fibre = [1 // _fiber_summand_weight(imV, i, V) for i in 1:rank(V)]
    w = vcat(w_from_base, w_from_fibre)
    
    nMarks = count(i -> i==v, dg.marks)
    vpEval = evaluate_vertex_polynomial(u, w, nMarks, dg.genus[v], VPs)
    res = res * vpEval
  end
  
  return res
end