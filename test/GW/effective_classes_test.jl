@testset "Effective classes for higher codimension" begin
    G = projective_space(GKMGraph, 2)
    H2 = GKM_second_homology(G)
    c1 = first_chern_class(G)
    line = curve_class(G, first(GKMtools.edges(G)))
    effective = GKMtools._effective_classes_with_functional_value

    # Divisor calls still select their original functional value.
    @test collect(effective(H2, c1, 3)) == [line]
    @test isempty(collect(effective(H2, c1, 1)))

    # Both endpoints of 0 <= codim(class) - c1*beta <= dim(G).
    @test collect(effective(H2, c1^2, 2)) == [zero(H2.H2)]
    @test collect(effective(H2, c1^3, 0)) == [line]
    @test collect(effective(H2, c1^5, 2)) == [line]
    @test isempty(collect(effective(H2, c1^3, 1)))
    @test isempty(collect(effective(H2, c1^3, -1)))
    @test isempty(collect(effective(H2, c1^3, 3)))
    @test collect(effective(H2, c1^2)) == [zero(H2.H2)]
    @test collect(effective(H2, c1^3)) == [line]
    @test collect(effective(H2, c1^6)) == [2 * line]
    @test_throws ArgumentError effective(H2, c1 + c1^2)
end
