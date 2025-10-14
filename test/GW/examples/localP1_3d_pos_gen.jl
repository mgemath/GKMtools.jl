# This code calculates the local GW invariants of rank 2 vector bundles over P1.

function local_p1_invariants_3d(a1, a2, dMax, gen)

  P1 = empty_gkm_graph(2, 3, ["p0", "p1"])
  g = gens(P1.M);
  add_edge!(P1, 1, 2, g[1]) 
  GMtoM = ModuleHomomorphism(P1.M, P1.M, [g[1], g[2], g[3]]);

  # The weights here don't matter as long as a1+a2 = -2. Otherwise they matter.
  w1 = g[2]
  w2 = -g[1] -g[2] # yields 1, -7, 55, ... for (-1, 3) (after *d^3)
  #w2 = g[2] - 2*g[1] # yields 1, 1, 1, ... for(-1, 3) (after *d^3). How is this possible? 
  #w2 = w2 + g[3]
  #w2 = g[3]

  w1p = w1 - a1*g[1]
  w2p = w2 - a2*g[1]

  V = vector_bundle(P1, P1.M, GMtoM, [w1 w2; w1p w2p])
  e = Edge(1, 2)
  V.con = Dict{Tuple{Edge, Int64}, Int64}((e, 1) => 1, (e, 2) => 2, (reverse(e), 1) => 1, (reverse(e), 2) => 2)
  b0 = curve_class(P1, Edge(1, 2))
  P_input = class_one()

  # Gromov-Witten invariants:
  res = Vector{Any}(undef, dMax)
  for d in 1:dMax
    res[d] = gromov_witten(V, d * b0, 0, P_input; show_bar=false, g=gen)
  end
  return res
end

###
### OUTPUT:
###

# julia> local_p1_invariants_3d(-1, -1, 4, 0, H)
# 4-element Vector{Any}:
#  1
#  1//8
#  1//27
#  1//64

# julia> local_p1_invariants_3d(-1, -1, 4, 1, H)
# 4-element Vector{Any}:
#  1//12
#  1//24
#  1//36
#  1//48

# julia> local_p1_invariants_3d(-1, -1, 4, 2, H)
# ^[[A^[[D4-element Vector{Any}:
#  1//240
#  1//120
#  1//80
#  1//60

# julia> local_p1_invariants_3d(-1, -1, 4, 3, H)
# 4-element Vector{Any}:
#  1//6048
#  1//756
#  1//224
#  2//189

# julia> local_p1_invariants_3d(-1, -1, 4, 4, H)
# 4-element Vector{Any}:
#  1//172800
#  1//5400
#  9//6400
#  4//675