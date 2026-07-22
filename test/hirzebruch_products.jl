@testset "Hirzebruch surface products" begin
  product_graph(G1, G2) = *(
    G1, G2;
    calculateCurveClasses=false,
    calculateConnection=false,
  )

  P1 = projective_space(GKM_graph, 1)
  systems = [
    system_of_invariants_6d(product_graph(
      gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, a)),
      P1,
    )) for a in 0:5
  ]

  comparison_modes = (
    (preserve_almost_complex=false,),
    (preserve_almost_complex=true,),
  )

  for a in 0:5, b in (a + 1):5, options in comparison_modes
    result = compare_systems(
      systems[a + 1], systems[b + 1];
      options...,
      primes=[2],
      integral_search_bounds=[1, 2, 3, 4],
    )
    @test result.status == (iseven(a - b) ? :equivalent : :not_equivalent)
    if result.status == :equivalent
      @test !isnothing(result.witness)
    else
      @test !isnothing(result.obstruction)
    end
  end
end
