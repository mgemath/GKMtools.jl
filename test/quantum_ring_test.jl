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
