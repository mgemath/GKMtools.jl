using Test, Oscar, GKMtools

function type_a_connection_profiles(G, v)
  C = get_connection(G)
  outgoing = [Edge(v, u) for u in all_neighbors(G.g, v)]
  return sort([Tuple(sort(Int.(C.a[e]))) for e in outgoing])
end

function all_type_a_connection_profiles(G)
  return sort([
    profile for v in vertices(G.g) for
    profile in type_a_connection_profiles(G, v)
  ])
end

@testset "Type A flag connections" begin
  @testset "Default type A flag connection is geometric" begin
    G_default = flag_variety(GKM_graph, [1, 1, 1])
    G_explicit = flag_variety(
      GKM_graph,
      [1, 1, 1];
      connection=:geometric,
    )

    C_default = get_connection(G_default)
    C_explicit = get_connection(G_explicit)
    @test C_default.a == C_explicit.a
    @test C_default.con == C_explicit.con
  end

  @testset "Connections on type A2 complete flags" begin
    G_geo = flag_variety(GKM_graph, [1, 1, 1]; connection=:geometric)
    G_comb = flag_variety(
      GKM_graph,
      [1, 1, 1];
      connection=:combinatorial,
    )

    @test isvalid(get_connection(G_geo))
    @test isvalid(get_connection(G_comb))

    for v in vertices(G_geo.g)
      @test type_a_connection_profiles(G_geo, v) == [
        (0, 0, 2),
        (0, 0, 2),
        (1, 1, 2),
      ]
      @test type_a_connection_profiles(G_comb, v) == [
        (-1, 1, 2),
        (-1, 1, 2),
        (1, 1, 2),
      ]
    end

    @test get_connection(G_geo).a != get_connection(G_comb).a
    @test get_connection(G_geo).con != get_connection(G_comb).con
  end

  @testset "Grassmannian connections agree" begin
    G_default = grassmannian(GKM_graph, 2, 4)
    G_geo = grassmannian(GKM_graph, 2, 4; connection=:geometric)
    G_comb = grassmannian(GKM_graph, 2, 4; connection=:combinatorial)

    C_default = get_connection(G_default)
    C_geo = get_connection(G_geo)
    C_comb = get_connection(G_comb)

    @test isvalid(C_default)
    @test isvalid(C_geo)
    @test isvalid(C_comb)
    @test C_default.a == C_geo.a == C_comb.a
    @test C_default.con == C_geo.con == C_comb.con

    for v in vertices(G_geo.g)
      @test type_a_connection_profiles(G_geo, v) == fill((0, 1, 1, 2), 4)
    end
  end

  @testset "Grassmannian connections agree through dimension 6" begin
    for b in 2:6
      for a in 1:(b - 1)
        G_geo = grassmannian(GKM_graph, a, b; connection=:geometric)
        G_comb = grassmannian(
          GKM_graph,
          a,
          b;
          connection=:combinatorial,
        )

        C_geo = get_connection(G_geo)
        C_comb = get_connection(G_comb)
        @test isvalid(C_geo)
        @test isvalid(C_comb)
        @test C_geo.a == C_comb.a
        @test C_geo.con == C_comb.con
      end
    end
  end

  @testset "Grassmannian forwards connection option" begin
    for option in (:geometric, :combinatorial)
      G_gr = grassmannian(GKM_graph, 2, 4; connection=option)
      G_fl = flag_variety(GKM_graph, [2, 2]; connection=option)

      @test get_connection(G_gr).a == get_connection(G_fl).a
      @test get_connection(G_gr).con == get_connection(G_fl).con
    end
  end

  @testset "Invalid type A connection option" begin
    @test_throws ArgumentError flag_variety(
      GKM_graph,
      [1, 1, 1];
      connection=:unknown,
    )
    @test_throws ArgumentError grassmannian(
      GKM_graph,
      2,
      4;
      connection=:unknown,
    )
  end

  @testset "Tangent directions and validity" begin
    examples = (([1, 1, 1], false), ([1, 2, 1], false), ([2, 2], true))
    for (s, is_grassmannian) in examples
      for option in (:geometric, :combinatorial)
        G = if is_grassmannian
          grassmannian(GKM_graph, 2, 4; connection=option)
        else
          flag_variety(GKM_graph, s; connection=option)
        end
        C = get_connection(G)
        @test isvalid(C)

        for base_edge in edges(G.g)
          for e in (base_edge, reverse(base_edge))
            source_index = G.edge_to_flag_index[e]
            target_index = G.edge_to_flag_index[reverse(e)]
            @test C.a[e][source_index] == ZZ(2)
            @test C.con[e][source_index] == target_index
          end
        end
      end
    end
  end

  @testset "Graph data are independent of connection option" begin
    G_geo = flag_variety(GKM_graph, [1, 1, 1]; connection=:geometric)
    G_comb = flag_variety(
      GKM_graph,
      [1, 1, 1];
      connection=:combinatorial,
    )

    @test collect(vertices(G_geo.g)) == collect(vertices(G_comb.g))
    @test collect(edges(G_geo.g)) == collect(edges(G_comb.g))
    @test G_geo.labels == G_comb.labels
    @test rank_torus(G_geo) == rank_torus(G_comb)
    @test base_ring(G_geo.M) == base_ring(G_comb.M)
    @test G_geo.weightType == G_comb.weightType

    for base_edge in edges(G_geo.g)
      for e in (base_edge, reverse(base_edge))
        @test GKMtools._w(G_geo, e) == GKMtools._w(G_comb, e)
      end
    end
  end

  @testset "Generalized flag profile cross-check" begin
    for option in (:geometric, :combinatorial)
      G_fast = flag_variety(GKM_graph, [1, 1, 1]; connection=option)
      G_general = generalized_gkm_flag(root_system(:A, 2); connection=option)
      @test all_type_a_connection_profiles(G_fast) ==
            all_type_a_connection_profiles(G_general)

      Gr_fast = grassmannian(GKM_graph, 2, 4; connection=option)
      Gr_general = generalized_gkm_flag(
        root_system(:A, 3),
        [1, 3];
        connection=option,
      )
      @test all_type_a_connection_profiles(Gr_fast) ==
            all_type_a_connection_profiles(Gr_general)
    end
  end
end
