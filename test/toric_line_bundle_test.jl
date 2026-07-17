using Test
using Oscar
using GKMtools

@testset "weighted projective line bundles" begin
  W = stacky_weighted_projective_space_fan([2, 3])
  L = line_bundle_O(W, 1)

  @test rank(L) == 1
  @test base_ring(GKMtools.lattice(baseof(L))) == QQ
  @test fiber_weight(L, 1, 1) == QQ(1, 3) * gens(L.M)[2]
  @test fiber_weight(L, 2, 1) == QQ(1, 2) * gens(L.M)[1]
  @test [Int(fiber_representation(L, v)[1, 1]) for v in vertices(baseof(L))] == [1, 1]

  Lminus = line_bundle_O(W, -1)
  @test [Int(fiber_representation(Lminus, v)[1, 1]) for v in vertices(baseof(Lminus))] == [2, 1]

  W1234 = stacky_weighted_projective_space_fan([1, 2, 3, 4])
  L2 = line_bundle_O(W1234, 2)
  @test any(!iszero(fiber_weight(L2, v, 1)) for v in vertices(baseof(L2)))
  @test !iszero(chern_class(L2, 1))
end
