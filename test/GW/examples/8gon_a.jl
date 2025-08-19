# Implementing the left hand GKM graph from [GKZ22, Example 2.44]
# As unsigned GKM graph, it is realized by the equivariant connected sum of three copies of S^2xS^2

G = gkm_2d([1 0; 0 1; -1 0; 0 -1; 1 0; 0 1; -1 0; 0 -1;])

H2 = GKM_second_homology(G)

S = QH_structure_constants(G)
pts = [QH_class(G, point_class(i, G)) for i in 1:8]
edg = [QH_class(G, poincare_dual(gkm_subgraph_from_vertices(G, [i, i%8 + 1]))) for i in 1:8]

# julia> QH_is_polynomial(G)
# QH is not polynomial for curve class (0, 1, 1, 0, 0, 0) at index CartesianIndex(4, 4, 4).
# false
# 
# julia> QH_is_homogeneous(G)
# true
# 
# julia> QH_is_commutative(G)
#true
#
#julia> QH_is_associative(G)
#Not associative for: (1, 1, 2)
#false

# (Q:) is this just because we didn't check high enough q coefficients?
# (A:) I guess calculating higher chern number stuff won't cancel out the non-associativity in the chern class 2 part!
# So this is actually NOT associative!

# Note: In any case a valid conclusion is that it is not Hamiltonian. But it could be that only the calculation of curve
# classes fails and associativity only fails because of a wrong assumption on relations in H_2.