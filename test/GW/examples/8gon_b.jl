# Implementing the right hand GKM graph from [GKZ22, Example 2.44]
# It does not admit a connection (or even an almost complex structure, I think), and cannot be realized by a Hamiltonian T-2 action.
# As unsigned GKM graph, it is realized by the equivariant connected sum of three copies of S^2xS^2

G = gkm_2d([1 0; 0 1; 1 0; 0 -1; -1 0; 0 1; 1 0; 0 -1])

H2 = GKM_second_homology(G) # crashes because c1 is not a GKM class. That's ok.

S = QH_structure_constants(G)