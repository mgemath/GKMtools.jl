# This code calculates the local GW invariants of rank 2 vector bundles over P1.

function local_p2(a, dMax, gen)

  P2 = empty_gkm_graph(3, 3, ["p0", "p1", "p2"])
  g = gens(P2.M);
  add_edge!(P2, 1, 2, g[1]) 
  add_edge!(P2, 2, 3, g[2])
  add_edge!(P2, 1, 3, g[1] + g[2])
  GMtoM = ModuleHomomorphism(P2.M, P2.M, [g[1]; g[2]; g[3]]);

  V = line_bundle(P2, P2.M, GMtoM, [g[3], -a*g[1]+g[3], -a*(g[1]+g[2])+g[3]])
  
  b0 = curve_class(P2, Edge(1, 2))
  P_input = class_one()

  # Gromov-Witten invariants:
  res = Vector{Any}(undef, dMax)
  for d in 1:dMax
    res[d] = gromov_witten(V, d * b0, 0, P_input; show_bar=false, g=gen)
  end
  return (V, res)
end