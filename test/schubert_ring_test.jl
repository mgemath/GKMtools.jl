@testset "Oscar graded Schubert cohomology ring" begin
  R = root_system(:A, 2)
  P2 = generalized_gkm_flag(R, [1])

  H, sigma = schubert_cohomology_ring(P2)
  h = only(sigma[2])
  point = only(sigma[3])
  @test is_graded(H)
  @test length.(sigma) == [1, 1, 1]
  @test vector_space_dimension(H) == 3
  @test h^2 == point
  @test h^3 == 0
  @test degree(h)[1] == 2
  @test degree(point)[1] == 4

  Hcodim, sigma_codim = schubert_cohomology_ring(
    P2;
    degree_convention=:codimension,
  )
  @test degree(only(sigma_codim[2]))[1] == 1
  @test vector_space_dimension(Hcodim) == 3

  flag = generalized_gkm_flag(R)
  Hflag, sigma_flag = schubert_cohomology_ring(flag)
  @test length.(sigma_flag) == [1, 2, 2, 1]
  @test vector_space_dimension(Hflag) == 6
  @test all(iszero(x^4) for x in sigma_flag[2])

  @test_throws ArgumentError schubert_cohomology_ring(P2; degree_convention=:unknown)
end
