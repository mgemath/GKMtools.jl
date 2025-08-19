# Implementing the GKM graph from [GKZ20, after Cor 2.13]
#
# It has a 2-torus acting and is 3-valent, so not 3-independent but still admits a unique connection.
# It is realized as smooth projective variety by P(E), where E -> X is an algebraic vector bundle over X = P^1 x P^1 (cf [GKZ20, Thm 5.1 and Prop 6.1])
# In terminology of [GKZ20, Def 4.9 & Cor 4.11], it has the maximum number of internal vertices.
# Hence, by [GKZ20, Thm 7.1 (I)], it cannot be realized by a complex structure with a compatible invariant Kaehler form (cf. [GKZ20, Remark 5.2]).
# So I guess the action on the algebraic variety P(E) is simply not algebraic...

g = Graph{Undirected}(8)
labels = ["v_$i" for i in 1:8]
M = free_module(ZZ, 2)
(t1, t2) = gens(M)
w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
G = gkm_graph(g, labels, M, w)

weights = [t1, t2, -t1, -t2]
for i in 1:8
  add_edge!(G, i, i%8 + 1, weights[(i-1) % 4 + 1])
end
add_edge!(G, 1, 5, t1+t2)
add_edge!(G, 2, 6, -t1+t2)
add_edge!(G, 3, 7, -t1-t2)
add_edge!(G, 4, 8, -t1-t2)

# julia> print_curve_classes(G)
# v_2 -> v_1: (0, 1, 1), Chern number: 4
# v_3 -> v_2: (1, 1, 0), Chern number: 4
# v_4 -> v_3: (0, 0, 1), Chern number: 2
# v_5 -> v_1: (0, 1, 0), Chern number: 2
# v_5 -> v_4: (1, 0, 0), Chern number: 2
# v_6 -> v_2: (0, 1, 0), Chern number: 2
# v_6 -> v_5: (0, -1, 1), Chern number: 0
# v_7 -> v_3: (0, 1, 0), Chern number: 2
# v_7 -> v_6: (1, -1, 0), Chern number: 0
# v_8 -> v_1: (1, 0, 0), Chern number: 2
# v_8 -> v_4: (0, 1, 0), Chern number: 2
# v_8 -> v_7: (0, 0, 1), Chern number: 2

b0a = curve_class(G, Edge(5, 6))
b0b = curve_class(G, Edge(6, 7))

println("Multiples of b0a = $b0a:")
for d in 1:5
  println("Class $d*b0a = $(d*b0a), 0 points -> $(gromov_witten(G, d*b0a, 0, class_one()))")
end
println("Multiples of b0b = $b0b:")
for d in 1:5
  println("Class $d*b0b = $(d*b0b), 0 points -> $(gromov_witten(G, d*b0b, 0, class_one()))")
end

# Multiples of b0a = (0, -1, 1):
# Class 1*b0a = (0, -1, 1), 0 points -> -1
# Class 2*b0a = (0, -2, 2), 0 points -> -1//8
# Class 3*b0a = (0, -3, 3), 0 points -> -1//27
# Class 4*b0a = (0, -4, 4), 0 points -> -1//64
# Class 5*b0a = (0, -5, 5), 0 points -> -1//125
# Multiples of b0b = (1, -1, 0):
# Class 1*b0b = (1, -1, 0), 0 points -> -1
# Class 2*b0b = (2, -2, 0), 0 points -> -1//8
# Class 3*b0b = (3, -3, 0), 0 points -> -1//27
# Class 4*b0b = (4, -4, 0), 0 points -> -1//64
# Class 5*b0b = (5, -5, 0), 0 points -> -1//125

# julia> for d in 1:5
#          println("Class $d*b0b = $(d*b0b), 0 points -> $(gromov_witten(G, d*b0b, 1, ev(1, point_class(6, G)))))")
#        end
# Class 1*b0b = (1, -1, 0), 0 points -> t1^2 - t1*t2)
# Class 2*b0b = (2, -2, 0), 0 points -> 1//4*t1^2 - 1//4*t1*t2)
# Class 3*b0b = (3, -3, 0), 0 points -> 1//9*t1^2 - 1//9*t1*t2)
# Class 4*b0b = (4, -4, 0), 0 points -> 1//16*t1^2 - 1//16*t1*t2)
 
# julia> for d in 1:5
#        println("Class $d*b0b = $(d*b0b), 0 points -> $(gromov_witten(G, d*b0b, 2, ev(2, point_class(7, G))*ev(1, point_class(6, G)))))")
#        end
# Class 1*b0b = (1, -1, 0), 0 points -> -t1^4 + t1^2*t2^2)
# Class 2*b0b = (2, -2, 0), 0 points -> -1//2*t1^4 + 1//2*t1^2*t2^2)
# Class 3*b0b = (3, -3, 0), 0 points -> -1//3*t1^4 + 1//3*t1^2*t2^2)
# Class 4*b0b = (4, -4, 0), 0 points -> -1//4*t1^4 + 1//4*t1^2*t2^2)

# julia> for d in 1:4
#        println("Class $d*b0b = $(d*b0b), 0 points -> $(gromov_witten(G, d*b0b, 4, ev(4, point_class(7, G))*ev(3, point_class(6, G))* ev(2, point_class(7, G))*ev(1, point_class(6, G)))))")
#        end
# Class 1*b0b = (1, -1, 0), 0 points -> -t1^8 + 2*t1^6*t2^2 - t1^4*t2^4)
# Class 2*b0b = (2, -2, 0), 0 points -> -2*t1^8 + 4*t1^6*t2^2 - 2*t1^4*t2^4)
# Class 3*b0b = (3, -3, 0), 0 points -> -3*t1^8 + 6*t1^6*t2^2 - 3*t1^4*t2^4)
# Class 4*b0b = (4, -4, 0), 0 points -> -4*t1^8 + 8*t1^6*t2^2 - 4*t1^4*t2^4)

# Generators (0, -1, 1), (0, 1, 0), (1, -1, 0) of Chern classes 0, 2, 0, resp.
H2gens = [curve_class(G, e) for e in [Edge(5, 6), Edge(3, 7), Edge(6, 7)]]

# This takes long, but one can just abort whenever, and all betas that have been calculated will be stored in G.QH_structure_constants.
for i in 0:3, j in 0:3, k in 0:3
  b = i*H2gens[1] + j*H2gens[2] + k*H2gens[3]
  QH_structure_constants(G, b)
end

# Result:
# We get polynomial QH structure constants in all curve classes computed so far.
# The nonzero ones calculated so far are:
#(0, -5, 5)
#(0, 1, 0)
#(1, -1, 0)
#(1, 0, 0)
#(1, 1, 0)
#(3, -3, 0)
#(0, -4, 4)
#(0, 0, 0)
#(2, -2, 0)
#
# Furthermore, the coefficients for beta=(0, -j, j) are constant for j in {2, 3, 4, 5}, so probably we get some q^beta/(1-q^beta) terms!
# The same holds for beta=(j, -j, 0) for j in {1, 2, 3}
# Note that the other classes can not appear infinitely often (if the structure constants are polynomial) because they have positive chern number.





#############################################
# OLD (MANUALLY BUILT CONNECTION):          #
#############################################

# H2 = GKM_second_homology(G)
# 
# conDict1 = Dict{Tuple{Edge, Edge}, Edge}(
#     (Edge(1, 2), Edge(1, 5)) => Edge(2, 6),
#     (Edge(1, 2), Edge(1, 8)) => Edge(2, 3),
#     (Edge(2, 3), Edge(2, 1)) => Edge(3, 4),
#     (Edge(2, 3), Edge(2, 6)) => Edge(3, 7),
#     (Edge(3, 4), Edge(3, 7)) => Edge(4, 8),
#     (Edge(3, 4), Edge(3, 2)) => Edge(4, 5),
#     (Edge(4, 5), Edge(4, 3)) => Edge(5, 6),
#     (Edge(4, 5), Edge(4, 8)) => Edge(5, 1),
#     (Edge(5, 6), Edge(5, 1)) => Edge(6, 2),
#     (Edge(5, 6), Edge(5, 4)) => Edge(6, 7),
#     (Edge(6, 7), Edge(6, 2)) => Edge(7, 3),
#     (Edge(6, 7), Edge(6, 5)) => Edge(7, 8),
#     (Edge(7, 8), Edge(7, 3)) => Edge(8, 4),
#     (Edge(7, 8), Edge(7, 6)) => Edge(8, 1),
#     (Edge(8, 1), Edge(8, 7)) => Edge(1, 2),
#     (Edge(8, 1), Edge(8, 4)) => Edge(1, 5),
#     # Now the diagonal edges:
#     (Edge(1, 5), Edge(1, 2)) => Edge(5, 6),
#     (Edge(1, 5), Edge(1, 8)) => Edge(5, 4),
#     (Edge(2, 6), Edge(2, 1)) => Edge(6, 5),
#     (Edge(2, 6), Edge(2, 3)) => Edge(6, 7),
#     (Edge(3, 7), Edge(3, 2)) => Edge(7, 6),
#     (Edge(3, 7), Edge(3, 4)) => Edge(7, 8),
#     (Edge(4, 8), Edge(4, 3)) => Edge(8, 7),
#     (Edge(4, 8), Edge(4, 5)) => Edge(8, 1)
# )
# for e in edges(G.g)
#   conDict1[(e, e)] = reverse(e)
# end
# for (e, ei) in keys(conDict1)
#   epi = conDict1[(e, ei)]
#   conDict1[reverse(e), epi] = ei
# end
# con1 = build_GKM_connection(G, conDict1)
# set_GKM_connection!(G, con1)