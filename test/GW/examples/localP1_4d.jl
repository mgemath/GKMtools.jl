# This code calculates the local GW invariants of rank 3 vector bundles over P1.

function local_p1_invariants_4d(a1, a2, a3, dMax)

  P1 = empty_gkm_graph(2, 4, ["p0", "p1"])
  g = gens(P1.M);
  add_edge!(P1, 1, 2, g[1]) 
  GMtoM = ModuleHomomorphism(P1.M, P1.M, [g[1], g[2], g[3], g[4]]);

  # The weights here don't matter as long as a1+a2 = -2. Otherwise they matter.
  w1 = g[2]
  w2 = g[3]
  w3 = -g[3] - g[2] - g[1]
  w3 = g[4]

  w1p = w1 - a1*g[1]
  w2p = w2 - a2*g[1]
  w3p = w3 - a3*g[1]

  V = vector_bundle(P1, P1.M, GMtoM, [w1 w2 w3; w1p w2p w3p])
  e = Edge(1, 2)
  V.con = Dict{Tuple{Edge, Int64}, Int64}((e, 1) => 1, (e, 2) => 2, (e, 3) => 3, (reverse(e), 1) => 1, (reverse(e), 2) => 2, (reverse(e), 3) => 3)
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

aMax = 6
dMax = 6
for a1 in -2:4, a2 in -2:4
  a3 = -2 - a1 - a2
  gw = local_p1_invariants_4d(a1, a2, a3, dMax)
  gw_normalized = [d^3 for d in 1:dMax] .* gw
  println("( GW(P1; $a1, $a2) )_d: ")
  for x in gw
    print(x)
    print(", ")
  end
  continue
  println("Times d^3: ")
  for x in gw_normalized
    print(x)
    print(", ")
  end
  n = a1 + 1
  if n != 0
    println("\n binomial($(n^2) d/ d)/$(n^2): ")
    println([binomial(n^2 * d, d) / n^2 for d in 1:dMax])
  else
    println("")
  end
  println()
end


