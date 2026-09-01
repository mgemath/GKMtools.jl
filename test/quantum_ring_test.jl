@testset "Oscar small quantum cohomology ring" begin
  P1 = projective_space(GKMGraph, 1)
  beta = curve_class(P1, Oscar.Edge(1, 2))
  QH, classes, q = small_quantum_cohomology_ring(P1; degrees=[beta])
  h = only(classes[2])
  q1 = only(q)

  @test is_graded(QH)
  @test length.(classes) == [1, 1]
  @test h^2 == q1
  @test degree(h)[1] == 2
  @test degree(q1)[1] == 4

  QHcodim, classes_codim, q_codim = small_quantum_cohomology_ring(
    P1; degrees=[beta], degree_convention=:codimension,
  )
  @test degree(only(classes_codim[2]))[1] == 1
  @test degree(only(q_codim))[1] == 2
  @test is_graded(QHcodim)

  @test_throws ArgumentError small_quantum_cohomology_ring(
    P1; degrees=[beta], degree_convention=:unknown,
  )
end

@testset "Threaded unmarked integration" begin
  P1 = projective_space(GKMGraph, 1)
  beta = 2 * curve_class(P1, Oscar.Edge(1, 2))
  products = [
    [point_class(P1, 1)],
    [first_chern_class(P1), point_class(P1, 2)],
  ]

  sequential = gromov_witten_nomarks(P1, beta, products; show_bar=false)
  threaded = gromov_witten_nomarks(
    P1, beta, products; show_bar=false, threaded=true,
  )
  @test threaded == sequential
end

@testset "Small quantum cohomology of the A2 flag variety" begin
  G = generalized_gkm_flag(root_system(:A, 2))
  beta = curve_class(G, Oscar.Edge(1, 2))

  QH, classes, q = small_quantum_cohomology_ring(G; degrees=[beta])

  @test is_graded(QH)
  @test length.(classes) == [1, 2, 2, 1]
  @test length(q) == 2

  sigma_s1, sigma_s2 = classes[2]
  sigma_s2s1, sigma_s1s2 = classes[3]
  @test sigma_s1^2 == sigma_s2s1 + q[1]
  @test sigma_s2^2 == sigma_s1s2
end
