using Test
using Oscar
using GKMtools

@testset "exceptional-divisor connections" begin
  affine_sp = affine_space(NormalToricVariety, 4)
  ambient = gkm_graph_of_toric(affine_sp)
  center = subgraph_from_vertices(ambient, [1])

  ordinary_blowup = blow_up(center)
  exceptional = subgraph(ordinary_blowup)
  @test !isempty(GKMtools.transport(GKMtools.connection(exceptional)))
  @test GKMtools.is_valid(exceptional, GKMtools.connection(exceptional))
  @test GKMtools.is_valid(
    ambient_graph(ordinary_blowup),
    GKMtools.connection(ambient_graph(ordinary_blowup)),
  )

  weighted_blowup = blow_up(center, [1, 2, 3, 4])
  weighted_exceptional = subgraph(weighted_blowup)
  @test base_ring(GKMtools.lattice(weighted_exceptional)) == QQ
  @test !isempty(GKMtools.transport(GKMtools.connection(weighted_exceptional)))
  @test GKMtools.is_valid(
    weighted_exceptional,
    GKMtools.connection(weighted_exceptional),
  )
  weighted_ambient = ambient_graph(weighted_blowup)
  @test !isempty(GKMtools.transport(GKMtools.connection(weighted_ambient)))
  @test GKMtools.is_valid(
    weighted_ambient,
    GKMtools.connection(weighted_ambient),
  )
end
