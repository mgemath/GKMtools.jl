G = GKMproj_space(4)
S1 = GKMsubgraph_from_vertices(G, ["x_0", "x_1", "x_2"]) # valid GKMsubgraph
S2 = GKMsubgraph_from_edges(G, [Edge(1, 2), Edge(2, 3), Edge(3,4), Edge(4, 1)]) # also valid GKM subgraph, but incompatible with connection C (see below)

PDS2 = PDClass(S2)
for i in 1:rank(G.equivariantCohomology.cohomRing)
  if PDS2[i] != 0
    println(factor(PDS2[i]))
  else
    println(0)
  end
end

@req pointClass(1, G) == PDClass(GKMsubgraph_from_vertices(G, [1])) "fail1"

@req GKM_isValid(G) "fail2"
@req GKM_isValidSubgraph(S1) "fail3"
@req GKM_isValidSubgraph(S2) "fail4"

C = get_GKM_connection(G)
@req GKM_isValidConnection(C) "fail5"
@req isCompatible(S1, C) "fail6"
@req !isCompatible(S2, C) "fail7"
