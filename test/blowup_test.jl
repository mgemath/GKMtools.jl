G = GKMproj_space(3)

S = GKMsubgraph_from_vertices(G, [1, 2])
blowupSub = blowupGKM(S)

Spoint = GKMsubgraph_from_vertices(G, [1])
blowupPt = blowupGKM(Spoint)