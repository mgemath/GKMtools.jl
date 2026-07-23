@testset "Tautological bundle from generalized GKM graph" begin
  R = root_system(:A, 3)
  lambda = fundamental_weight(R, 1)

  G = generalized_gkm_flag(R)
  E = tautological_bd(lambda, G)
  legacy_E = tautological_bd(lambda)
  @test rank(E) == rank(legacy_E)
  @test [fiber_weight(E, v, i) for v in vertices(baseof(E)), i in 1:rank(E)] ==
    [fiber_weight(legacy_E, v, i) for v in vertices(baseof(legacy_E)), i in 1:rank(legacy_E)]

  P = generalized_gkm_flag(R, [1, 3])
  bundles = tautological_bd([lambda, 2 * lambda], P)
  @test rank.(bundles) == rank.(tautological_bd([lambda, 2 * lambda], [1, 3]))

  R2 = root_system(:A, 2)
  @test_throws ArgumentError tautological_bd(
    fundamental_weight(R2, 1), G
  )
end
