using Test, Combinatorics, Oscar, GKMtools

@testset "GKM costructions" begin

  M = free_module(ZZ, 2);
  g = Graph{Undirected}(3);
  add_edge!(g, 1, 2);
  add_edge!(g, 2, 3);
  add_edge!(g, 1, 3);
  de = Dict(edges(g) .=> [zero(M) for _ in n_edges(g)]);
  labs = ["a", "b", "c"];
  G = gkm_graph(g, labs, M, de) 

  @test rank_torus(G) == 2
  @test valency(G) == 2
  @test typeof(convert_weights(G)) == GKM_graph{QQFieldElem}

#   @testset "Betti numbers of Hirzebruch surfaces" begin

#   end


end