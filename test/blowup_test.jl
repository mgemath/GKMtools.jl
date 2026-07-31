using Test

G = GKMproj_space(3)

S = GKMsubgraph_from_vertices(G, [1, 2])
blowupSub = blowupGKM(S)

Spoint = GKMsubgraph_from_vertices(G, [1])
blowupPt = blowupGKM(Spoint)

@test other_vertex(G, 1, 1) == 2
@test other_vertex(G, 1, 2) == 3
@test_throws ArgumentError other_vertex(G, 1, 4)

affine_sp = affine_space(NormalToricVariety, 4)
affine_graph = gkm_graph_of_toric(affine_sp)
affine_subgraph = subgraph_from_vertices(affine_graph, [1])
@test_nowarn blow_up(affine_subgraph, [1, 2, 3, 4])