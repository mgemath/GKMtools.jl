G = GKMproj_space(3)
t1, t2, t3, t4 = gens(G.equivariantCohomology.coeffRing)
H = [t1, t2, t3, t4]
p = QH_class(G, H)
println((p-t1)*(p-t2)*(p-t3)*(p-t4)) # should be q^(1)

S = GKMsubgraph_from_vertices(G, [1, 2])
s = PDClass(S)