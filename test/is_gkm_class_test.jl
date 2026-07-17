using Test
using Oscar
using GKMtools

@testset "GKM class divisibility" begin
  G = projective_space(GKMGraph, 1)
  e = gens_cohomRing(G)
  t = gens_coeffRing(G)

  @test is_gkm_class(G, one(e[1]))
  @test is_gkm_class(G, weight_class(G, first(edges(G))) * e[1])

  # Localizations (1, 0) do not differ by a multiple of the edge weight.
  @test !is_gkm_class(G, e[1])

  # Equal rational localizations satisfy the edge congruence only after
  # localization, but do not belong to the polynomial GKM ring.
  @test !is_gkm_class(G, inv(t[1]) * one(e[1]))
end
