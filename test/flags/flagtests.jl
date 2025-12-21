test_blowup = true # keep this for later.

G = empty_gkm_graph(3, 2, ["1", "2", "3"])
g1, g2 = gens(G.M)
add_edge!(G, 1, 2, g1)
add_edge!(G, 2, 3, -g2)
add_standalone_flag!(G, 1, g2)
add_standalone_flag!(G, 1, -g1-g2)
add_standalone_flag!(G, 2, g1+g2)
add_standalone_flag!(G, 3, g1)
add_standalone_flag!(G, 3, -g1-g2)

@req isvalid(G) "G is invalid"

C = get_any_connection(G)
@req isvalid(C) "Connection for G is invalid"

set_connection!(G, get_any_connection(G))
S_cpct = gkm_subgraph_from_vertices(G, [1, 2])
S_non_cpct = gkm_subgraph_from_vertices(G, [1, 2]; include_standalone_flags=true)
@req isvalid(S_cpct) "S_cpct is invalid"
@req isvalid(S_non_cpct) "S_non_cpct is invalid"

if test_blowup
  # Problem with blowup: Doesn't produce the standalone flags.
  B1 = blow_up(S_cpct)
  @req isvalid(B1) "Blowup of S_cpct is invalid"
end

# Test products:

P1 = projective_space(GKM_graph, 1)
P = G * P1
@req isvalid(P) "product is invalid"

# Test subspace and blowup again:
S = gkm_subgraph_from_vertices(P, ["1,1", "1,2", "2,2"]; include_standalone_flags=true)
@req isvalid(S) "S is not valid"
if test_blowup
  try
    BS = blow_up(S)
    @req isvalid(BS) "Blowup of S is invalid"
  catch e
    if isa(e, ArgumentError) && occursin("constant codimension", e.msg)
      println("Skipping blowup of S: does not have constant codimension")
    else
      rethrow(e)
    end
  end
end

# Test total space and its connection properties
F = gkm_3d_twisted_flag()
set_connection!(F, get_any_connection(F))
TF = tangent_bd(F)
TF_tot = total_space(TF)

set_connection!(G, get_any_connection(G))
TG = tangent_bd(G)
TG_tot = total_space(TG)