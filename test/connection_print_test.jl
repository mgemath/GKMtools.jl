using Test, Oscar, GKMtools

@testset "Stable connection printing" begin
  G = generalized_gkm_flag(root_system(:A, 2))
  output = sprint(io -> print_connection(io, G))
  for section in split(output, r"Connection:\n|a_i's:\n")[2:end]
    rows = filter(!isempty, split(section, '\n'))
    endpoints = map(rows) do row
      source, target = split(first(split(row, " => ")), " -> ")
      (number_vertex_of_label(G, String(source)), number_vertex_of_label(G, String(target)))
    end
    @test endpoints == sort(endpoints)
  end
  # Rebuild both dictionaries in the opposite insertion order.
  for data in (G.connection.transport, G.connection.a)
    entries = reverse(collect(data))
    empty!(data)
    for (edge, value) in entries
      data[edge] = value
    end
  end
  @test sprint(io -> print_connection(io, G)) == output
  @test sprint(io -> print_connection(io, G; verbose=false)) ==
    "\nBirkhoff-Grothendieck connection for GKM graph with 6 nodes and valency 3"
end
