# This code calculates the local GW invariants of rank 2 vector bundles over P1.

function local_p1_invariants_3d(a1, a2, dMax, gen, H)

  P1 = empty_gkm_graph(2, 3, ["p0", "p1"])
  g = gens(P1.M);
  add_edge!(P1, 1, 2, g[1]) 
  GMtoM = ModuleHomomorphism(P1.M, P1.M, [g[1], g[2], g[3]]);

  # The weights here don't matter as long as a1+a2 = -2. Otherwise they matter.
  w1 = g[2]
  w2 = -g[1] -g[2] # yields 1, -7, 55, ... for (-1, 3) (after *d^3)
  #w2 = g[2] - 2*g[1] # yields 1, 1, 1, ... for(-1, 3) (after *d^3). How is this possible? 
  #w2 = w2 + g[3]

  w1p = w1 - a1*g[1]
  w2p = w2 - a2*g[1]

  V = vector_bundle(P1, P1.M, GMtoM, [w1 w2; w1p w2p])
  e = Edge(1, 2)
  V.con = Dict{Tuple{Edge, Int64}, Int64}((e, 1) => 1, (e, 2) => 2, (reverse(e), 1) => 1, (reverse(e), 2) => 2)
  b0 = curve_class(P1, Edge(1, 2))
  P_input = class_one()

  # Gromov-Witten invariants:
  res = Vector{Any}(undef, dMax)
  contribs = Vector{Any}(undef, dMax)
  coeffs = Vector{Any}(undef, dMax)
  for d in 1:dMax
    println()
    println("Degree $d:")
    res[d], contribs[d] = gromov_witten_pos_gen(V, d * b0, 0, gen, [P_input], H; show_bar=false)
    # println("Contribs: $(contribs[d])")
    # _find_coefficients(contribs[d], d)
  end
  return res
end

function _find_coefficients(l::Vector, d::Int64)
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  l = l .* lcm(denominator.(l))
  l = [evaluate(e, [x, one(x), zero(x)]) for e in l]
  l = [numerator(e) * (1 // coeff(denominator(e), 1)) for e in l] # should be QQMPolyRingElem now.
  M = matrix(QQ, [coeff(e, x^i) for e in l, i in 0:2*d])
  K = kernel(M; side=:left)
  
  println("M=$M\nK=$K")
  s = solve(M, vcat(zeros(QQ, 2*d), [QQ(1)//(12*d)]))
  println("s=$s")
end

