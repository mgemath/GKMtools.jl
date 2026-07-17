using Test, Oscar, GKMtools

function connection_profiles_at_vertex(G, C, v)
  outgoing = [Edge(v, u) for u in all_neighbors(G.g, v)]
  return sort([Tuple(sort(Int.(C.a[e]))) for e in outgoing])
end

@testset "Generalized flag connections" begin
  R = root_system(:A, 2)

  @testset "Default G/P connection is geometric" begin
    G_default = generalized_gkm_flag(R)
    G_explicit = generalized_gkm_flag(R; connection=:geometric)

    C_default = get_connection(G_default)
    C_explicit = get_connection(G_explicit)

    @test C_default !== nothing
    @test C_explicit !== nothing
    @test C_default.a == C_explicit.a
    @test C_default.con == C_explicit.con
  end

  @testset "Geometric connection on SL(3)/B" begin
    G = generalized_gkm_flag(R; connection=:geometric)
    C = get_connection(G)

    @test C !== nothing
    @test isvalid(C)

    id_vertex = findfirst(==("id"), G.labels)
    @test id_vertex !== nothing
    @test connection_profiles_at_vertex(G, C, id_vertex) == [
      (0, 0, 2),
      (0, 0, 2),
      (1, 1, 2),
    ]
  end

  @testset "Combinatorial connection on SL(3)/B" begin
    G = generalized_gkm_flag(R; connection=:combinatorial)
    C = get_connection(G)

    @test C !== nothing
    @test isvalid(C)

    id_vertex = findfirst(==("id"), G.labels)
    @test id_vertex !== nothing
    @test connection_profiles_at_vertex(G, C, id_vertex) == [
      (-1, 1, 2),
      (-1, 1, 2),
      (1, 1, 2),
    ]
  end

  @testset "Invalid G/P connection option" begin
    @test_throws ArgumentError generalized_gkm_flag(R; connection=:unknown)
  end

  @testset "Diagonal entries and edge transport" begin
    for connection in (:geometric, :combinatorial)
      G = generalized_gkm_flag(R; connection=connection)
      C = get_connection(G)

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

  @testset "Parabolic overloads" begin
    for parabolic_indices in ([1], [2])
      G_geo = generalized_gkm_flag(
        R,
        parabolic_indices;
        connection=:geometric,
      )
      C_geo = get_connection(G_geo)

      @test isvalid(C_geo)
      for v in vertices(G_geo.g)
        profiles = connection_profiles_at_vertex(G_geo, C_geo, v)
        @test all(profile == (1, 2) for profile in profiles)
      end

      G_comb = generalized_gkm_flag(
        R,
        parabolic_indices;
        connection=:combinatorial,
      )
      @test isvalid(get_connection(G_comb))
    end

    simple_parabolic = [simple_root(R, 1)]
    for connection in (:geometric, :combinatorial)
      G = generalized_gkm_flag(R, simple_parabolic; connection=connection)
      @test isvalid(get_connection(G))
    end
  end

  @testset "Root-system consistency" begin
    examples = [
      (root_system(:B, 2), Int[]),
      (root_system(:B, 2), [1]),
      (root_system(:G, 2), Int[]),
      (root_system(:G, 2), [1]),
    ]

    for (root_system_example, parabolic_indices) in examples
      for connection in (:geometric, :combinatorial)
        G = generalized_gkm_flag(
          root_system_example,
          parabolic_indices;
          connection=connection,
        )
        @test isvalid(get_connection(G))
      end

      G_geo = generalized_gkm_flag(
        root_system_example,
        parabolic_indices;
        connection=:geometric,
      )
      @test all(
        aij >= 0 for coefficients in values(get_connection(G_geo).a) for aij in coefficients
      )
    end
  end
end
