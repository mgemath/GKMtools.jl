using Test
using Oscar
using GKMtools

@testset "weighted projective line bundles" begin
  W = stacky_weighted_projective_space_fan([2, 3])
  @test gkm_graph_of_orbifold_toric(W, small_torus=false) == gkm_graph_of_orbifold_toric(W, small_torus=false)
  @test gkm_graph_of_orbifold_toric(W, small_torus=false) != gkm_graph_of_orbifold_toric(stacky_weighted_projective_space_fan([1, 2]), small_torus=false)
  L = line_bundle_O(W, 1)

  @test rank(L) == 1
  @test base_ring(GKMtools.lattice(baseof(L))) == QQ
  @test fiber_weight(L, 1, 1) == -QQ(1, 3) * gens(L.M)[2]
  @test fiber_weight(L, 2, 1) == -QQ(1, 2) * gens(L.M)[1]
  @test [Int(fiber_representation(L, v)[1, 1]) for v in vertices(baseof(L))] == [1, 1]

  Lr2 = line_bundle_O(W, [1, 0])
  Lr1 = line_bundle_O(W, [0, 1])
  @test [fiber_weight(Lr2, v, 1) for v in vertices(baseof(Lr2))] == [fiber_weight(line_bundle_O(W, 2), v, 1) for v in vertices(baseof(Lr2))]
  @test [fiber_weight(Lr1, v, 1) for v in vertices(baseof(Lr1))] == [fiber_weight(line_bundle_O(W, 3), v, 1) for v in vertices(baseof(Lr1))]
  @test Lr2 == line_bundle_O(W, 2)
  @test Lr1 == line_bundle_O(W, 3)
  @test Lr2 != Lr1

  Lminus = line_bundle_O(W, -1)
  @test [Int(fiber_representation(Lminus, v)[1, 1]) for v in vertices(baseof(Lminus))] == [2, 1]

  W1234 = stacky_weighted_projective_space_fan([1, 2, 3, 4])
  Ldiv = line_bundle_O(W1234, [1, 0, 1, -1])
  Ldegree = line_bundle_O(W1234, 1 + 3 - 4)
  @test [fiber_weight(Ldiv, v, 1) for v in vertices(baseof(Ldiv))] == [fiber_weight(Ldegree, v, 1) for v in vertices(baseof(Ldegree))]

  L2 = line_bundle_O(W1234, 2)
  @test any(!iszero(fiber_weight(L2, v, 1)) for v in vertices(baseof(L2)))
  @test !iszero(chern_class(L2, 1))

  @test integrate(baseof(L), chern_class(L, 1)) == 1//6

  Lr1 = line_bundle_O(W, [0, 1]) # O(3)
  @test integrate(baseof(Lr1), chern_class(Lr1, 1)) == 1//2

  Lr2 = line_bundle_O(W, [1, 0]) # O(2)
  @test integrate(baseof(Lr2), chern_class(Lr2, 1)) == 1//3

  Lminus = line_bundle_O(W, -1)
  @test integrate(baseof(Lminus), chern_class(Lminus, 1)) == -1//6
end
