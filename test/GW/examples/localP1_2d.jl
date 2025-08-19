# This code calculates the local GW invariants of line bundles over P1.

function local_p1_invariants_2d(a, dMax)

  P1 = empty_gkm_graph(2, 2, ["p0", "p1"])
  g = gens(P1.M);
  add_edge!(P1, 1, 2, g[1]) 
  GMtoM = ModuleHomomorphism(P1.M, P1.M, [g[1], g[2]]);

  w1 = g[2]

  w1p = w1 - a*g[1]

  V = line_bundle(P1, P1.M, GMtoM, [w1, w1p])
  e = Edge(1, 2)
  V.con = Dict{Tuple{Edge, Int64}, Int64}((e, 1) => 1, (reverse(e), 1) => 1)
  b0 = curve_class(P1, Edge(1, 2))
  P_input = class_one()

  # Gromov-Witten invariants:
  res = Vector{Any}(undef, dMax)
  for d in 1:dMax
    res[d] = gromov_witten(V, d * b0, 0, P_input; show_bar=false)
  end
  return res
end

# The upshot of the following experiment is:
# Conjecture:
# GW_d(a, -2-a) = 1/d^3 * binomial(k^2 n, n) / k^2 up to sign, where k = |a+1|

dMax = 6
for a in -2:4
  gw = local_p1_invariants_2d(a, dMax)
  gw_normalized = [d^3 for d in 1:dMax] .* gw
  println("( GW(P1; $a) )_d: ")
  # for x in gw
  #   print(x)
  #   print(", ")
  # end
  println("Times d^3: ")
  for x in gw_normalized
    print(x)
    print(", ")
  end
  println()
end


