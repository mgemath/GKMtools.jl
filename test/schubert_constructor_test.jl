using Test, Oscar, GKMtools

@testset "Schubert variety smoothness diagnostics" begin
  @testset "Non-rationally smooth varieties" begin
    R = root_system(:A, 3)
    G = generalized_gkm_flag(R)
    singular_label = "s2*s1*s3*s2"
    singular_index = number_vertex_of_label(G, singular_label)

    # Test all public entry points, including the opposite Schubert interval.
    constructors = [
      () -> generalized_gkm_schubert(G, singular_label),
      () -> generalized_gkm_schubert(G, singular_index),
      () -> generalized_gkm_schubert(R, singular_label),
      () -> generalized_gkm_schubert(R, Int[], singular_label),
      () -> generalized_gkm_schubert(R, RootSpaceElem[], singular_label),
      () -> generalized_gkm_schubert(G, "s1*s3"; descending=false),
      () -> generalized_gkm_schubert(R, "s1*s3", false),
    ]
    for construct in constructors
      err = try
        construct()
      catch err
        err
      end
      @test err isa ArgumentError
      message = sprint(showerror, err)
      @test occursin("not rationally smooth", message)
      @test occursin("does not have constant valency", message)
      @test occursin("supports smooth Schubert varieties only", message)
    end
  end

  @testset "Rationally smooth but singular varieties" begin
    G = generalized_gkm_flag(root_system(:B, 2))
    err = try
      generalized_gkm_schubert(G, "s2*s1*s2")
    catch err
      err
    end
    @test err isa ArgumentError
    message = sprint(showerror, err)
    @test occursin("s2*s1*s2", message)
    @test occursin("rational smoothness does not imply smoothness", message)
    @test !occursin("not rationally smooth", message)
    @test occursin("Underlying error:", message)
    @test occursin("Chern weights are not proportional to the edge weight", message)
  end

  @testset "Smooth varieties" begin
    R = root_system(:B, 2)
    G = generalized_gkm_flag(R)
    S = generalized_gkm_schubert(G, "s1*s2")
    @test ambient_graph(S) === G
    @test GKMtools.num_vertices(subgraph(S)) == 4
    @test valency(subgraph(S)) == 2
    @test GKMtools.is_valid(subgraph(S), GKMtools.connection(subgraph(S)))

    # Smooth points, a partial flag variety, and opposite Schubert varieties.
    point = generalized_gkm_schubert(G, "id")
    @test GKMtools.num_vertices(subgraph(point)) == 1
    @test valency(subgraph(point)) == 0
    partial = generalized_gkm_schubert(R, [1], "s2")
    @test GKMtools.num_vertices(subgraph(partial)) == 2
    opposite = generalized_gkm_schubert(G, "s1*s2"; descending=false)
    @test GKMtools.num_vertices(subgraph(opposite)) == 4
    @test valency(subgraph(opposite)) == 2
  end
end
