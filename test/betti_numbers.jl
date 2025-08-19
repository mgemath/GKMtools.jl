@testset "GKM graphs-Betti numbers" begin

  @testset "Betti numbers of Hirzebruch surfaces" begin
    H6 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 6))
    @test betti_numbers(H6) == ZZ.([1, 2, 1])
    @test betti_number(H6, 0) == ZZ(1)
    @test betti_number(H6, 1) == ZZ(0)
    @test betti_number(H6, 2) == ZZ(2)
  end

  @testset "Betti numbers of blow-ups" begin
    B = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 6))
    @test betti_numbers(H6) == ZZ.([1, 2, 1])
    @test betti_number(H6, 0) == ZZ(1)
    @test betti_number(H6, 1) == ZZ(0)
    @test betti_number(H6, 2) == ZZ(2)
  end


end