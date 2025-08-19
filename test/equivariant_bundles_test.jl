G = GKMproj_space(2)
M = free_module(ZZ, 4)
g = gens(M)
GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]])
w1 = [g[1], g[2], g[3]]
V1 = line_bundle(G, M, GMtoM, w1)
w2 = [g[4], g[4], g[4]]
V2 = line_bundle(G, M, GMtoM, w2)
V = direct_sum(V1, V2)

B = projective_bundle(V)