using Test, Oscar, GKMtools
import Serialization

_gkm_coordinates(v) = [v[i] for i in 1:rank(parent(v))]

@testset "Generalized flag graph persistence" begin
  mktempdir() do dir
    path = joinpath(dir, "graph.jls")
    for (family, r, levi, method) in [
      (:A, 2, [1], :birkhoff_grothendieck),
      (:B, 2, Int[], :cartan),
      (:F, 4, [1, 2, 3], :cartan),
      (:A, 1, [1], :algorithm),
    ]
      G = generalized_gkm_flag(root_system(family, r), levi; connection=method)
      @test serialize_gkm_graph(path, G) === G
      d = Serialization.deserialize(path)
      @test d.version == 1
      H = deserialize_gkm_graph(path)
      @test typeof(H) == typeof(G)
      @test GKMtools.get_string.(H.core.labels) == GKMtools.get_string.(G.core.labels)
      @test [word(GKMtools.flag(v)) for v in H.core.labels] ==
        [word(GKMtools.flag(v)) for v in G.core.labels]
      @test H.core.edge_flags == G.core.edge_flags
      @test [[_gkm_coordinates(f.weight) for f in fs] for fs in H.core.flags] ==
        [[_gkm_coordinates(f.weight) for f in fs] for fs in G.core.flags]
      @test H.connection.transport == G.connection.transport
      @test H.connection.a == G.connection.a
      @test H.connection.connection_type == G.connection.connection_type
      @test H.H2.edge_to_gen == G.H2.edge_to_gen
      @test [_gkm_coordinates(H.H2.quotient(g)) for g in gens(H.H2.edge_lattice)] ==
        [_gkm_coordinates(G.H2.quotient(g)) for g in gens(G.H2.edge_lattice)]
      @test [_gkm_coordinates(H.H2.chern(g)) for g in gens(H.H2.H2)] ==
        [_gkm_coordinates(G.H2.chern(g)) for g in gens(G.H2.H2)]
      @test H.H2.ray_sum == G.H2.ray_sum
      if family == :A && r == 2
        @test quantum_schubert_table(quantum_schubert_context(H)) ==
          quantum_schubert_table(quantum_schubert_context(G))
      end
    end
    Serialization.serialize(path, (format="GKMtools.generalized_flag", version=2))
    @test_throws ArgumentError deserialize_gkm_graph(path)
    Serialization.serialize(path, [1, 2])
    @test_throws ArgumentError deserialize_gkm_graph(path)
  end
end
