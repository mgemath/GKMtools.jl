using Test
using Oscar
using GKMtools

@testset "stacky cones and fans" begin
  cone = StackyCone([2 0; 0 3])
  @test dim(cone) == 2
  @test n_rays(cone) == 2
  @test is_orbifold(cone)
  @test cone.primitive_rays_matrix == matrix(ZZ, [1 0; 0 1])
  @test_throws ArgumentError StackyCone([1 0; 2 0])
  @test_throws ArgumentError StackyCone([0 0])
  @test_throws ArgumentError StackyCone([1 0; 0 1]; ambient_rank=3)

  fan = StackyFan([2 0; 0 3; -1 -1], [[1, 2], [2, 3], [3, 1]])
  @test dim(fan) == 2
  @test n_rays(fan) == 3
  @test is_orbifold(fan)
  @test length(maximal_cones(fan)) == 3
  @test maximal_cones(fan)[2].ray_indices == [2, 3]
  @test sum(GKMtools.incidence_matrix(fan), dims=2) == fill(2, 3, 1)
  graph = gkm_graph_of_orbifold_toric(fan)
  @test nv(graph.core.g) == 3
  @test ne(graph.core.g) == 3

  ineffective = StackyFan([1 0; 0 1; -1 -1], [[1, 2], [2, 3], [3, 1]]; generic_stabilizer=[2])
  @test ineffective.generic_stabilizer == [2]
  ineffective_graph = gkm_graph_of_orbifold_toric(ineffective)
  @test all(I -> I.isotropy_group == [2], ineffective_graph.vertex_isotropy)
  @test all(I -> I.isotropy_group == [2], Iterators.flatten(ineffective_graph.flag_isotropy))

  fan_wps = stacky_weighted_projective_space_fan([1, 2, 3, 4])
  graph_wps = gkm_graph_of_orbifold_toric(fan_wps)
  @test nv(graph_wps.core.g) == 4

  @test_throws ArgumentError StackyFan([1 0; 0 1], [[1]])
  @test_throws ArgumentError StackyFan([1 0; 2 0], [[1, 2]])
end
