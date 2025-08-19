# This point calculates the famous numbers N_d of degree d rational curves in P^2 through 3d-1 points

P2 = projective_space(GKM_graph, 2)
p = point_class(P2, 2)
P_input = class_one()
beta = curve_class(P2, Edge(1, 2))

for d in 1:4 #1, 2, 3 are reasonably quick (and correct: 1, 1, 12), 4 takes lang already (should be ).
  if d == 1
    global P_input = P_input * ev(1, p) * ev(2, p)
  else
    global P_input = P_input * ev(3*d-3, p) * ev(3*d-2, p) * ev(3*d-1, p)
  end
  g = gromov_witten(P2, d*beta, 3*d-1, P_input; show_bar=true)
  println("N_d for d=$d -> $g")
end