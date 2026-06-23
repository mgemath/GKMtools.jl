using Oscar, GKMtools

# ─── helpers ──────────────────────────────────────────────────────────────────

# The scalar λ with pos(dst(e)) - pos(src(e)) = λ · w_e for a valid drawing.
function edge_lambda(G, positions, e)
  d = rank_torus(G)
  we = [QQ(GKMtools._w(G, e)[k]) for k in 1:d]
  delta = positions[dst(e)] - positions[src(e)]
  kstar = findfirst(k -> we[k] != 0, 1:d)
  return delta[kstar] // we[kstar]
end

# Count chambers among admissible_drawing_representatives whose per-edge sign
# vector is all-positive or (globally flipped) all-negative.
function count_positive_chambers(G)
  reps = admissible_drawing_representatives(G)
  cnt = 0
  for pos in reps
    lams = [edge_lambda(G, pos, e) for e in edges(G.g)]
    if all(>(0), lams) || all(<(0), lams)
      cnt += 1
    end
  end
  return cnt
end

# Hexagon 1-2-3-4-5-6-1 with weights x,y,x,y,x,y. Going around the cycle the
# x-coefficients must sum to zero, which is impossible with all multipliers
# positive, so this graph has no positive drawing.
function hexagon_xyxyxy()
  G = empty_gkm_graph(6, 2, ["$i" for i in 1:6])
  x, y = gens(G.M)
  add_edge!(G, 1, 2, x); add_edge!(G, 2, 3, y); add_edge!(G, 3, 4, x)
  add_edge!(G, 4, 5, y); add_edge!(G, 5, 6, x); add_edge!(G, 6, 1, y)
  return G
end

# ─── Test 1: P^2 has a positive drawing ───────────────────────────────────────
println("=== Test 1: P^2 = flag_variety([2,1]) ===")
G = flag_variety(GKM_graph, [2, 1])
res = positive_drawing_representative(G)
@assert res !== nothing "Expected a positive drawing for P^2"
@assert res.admissible "Positive representative must be admissible"
@assert length(res.positions) == n_vertices(G.g)
for e in edges(G.g)
  λ = edge_lambda(G, res.positions, e)
  @assert λ > 0 "Edge $e is not a positive multiple of its weight (λ = $λ)"
end
println("P^2 positive drawing found; convex = $(res.convex).")

# ─── Test 2: no positive drawing (empty positive chamber) ─────────────────────
println("\n=== Test 2: hexagon x,y,x,y,x,y (no positive drawing) ===")
Gh = hexagon_xyxyxy()
@assert positive_drawing_representative(Gh) === nothing "Expected NO positive drawing for the hexagon"
println("Hexagon: positive_drawing_representative returned nothing, as expected.")

# ─── Test 3: consistency with the chamber enumeration ─────────────────────────
# Validates the conjecture "at most one positive representative" and that
# existence agrees with the all-positive chamber actually being realizable.
println("\n=== Test 3: consistency with admissible_drawing_representatives ===")
for (name, H) in [("P^2", flag_variety(GKM_graph, [2, 1])),
                  ("P^3", flag_variety(GKM_graph, [1, 3])),
                  ("F(1,1,1)", flag_variety(GKM_graph, [1, 1, 1])),
                  ("hexagon", hexagon_xyxyxy())]
  cnt = count_positive_chambers(H)
  exists = positive_drawing_representative(H) !== nothing
  @assert cnt <= 1 "At most one positive chamber expected for $name, got $cnt"
  @assert (cnt == 1) == exists "Existence mismatch for $name (cnt=$cnt, exists=$exists)"
  println("  $name: positive chambers = $cnt, positive_drawing_representative exists = $exists")
end

# ─── Test 4: strong kwarg ─────────────────────────────────────────────────────
println("\n=== Test 4: strong kwarg ===")
rw = positive_drawing_representative(G; strong=false)
rs = positive_drawing_representative(G; strong=true)
@assert rw.convex isa Bool && rs.convex isa Bool
@assert rw.positions == rs.positions "strong kwarg must not change the representative"
println("P^2: weakly_convex = $(rw.convex), strongly_convex = $(rs.convex)")

# ─── Test 5: Bennis example (no convex drawing) — smoke test ───────────────────
println("\n=== Test 5: Bennis example smoke test ===")
Gb = empty_gkm_graph(12, 2, ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"])
let (x, y) = gens(Gb.M)
  add_edge!(Gb, "A", "B", x);   add_edge!(Gb, "B", "C", x-y); add_edge!(Gb, "C", "D", -y)
  add_edge!(Gb, "D", "E", -x);  add_edge!(Gb, "E", "F", y-x); add_edge!(Gb, "F", "A", y)
  add_edge!(Gb, "G", "H", x);   add_edge!(Gb, "H", "I", x-y); add_edge!(Gb, "I", "J", -y)
  add_edge!(Gb, "J", "K", -x);  add_edge!(Gb, "K", "L", y-x); add_edge!(Gb, "L", "G", y)
  add_edge!(Gb, "B", "H", x+y); add_edge!(Gb, "C", "I", x+y); add_edge!(Gb, "E", "K", x+y)
  add_edge!(Gb, "F", "L", x+y); add_edge!(Gb, "A", "J", x+y); add_edge!(Gb, "D", "G", x+y)
end
resb = positive_drawing_representative(Gb)
if resb === nothing
  println("Bennis: no positive drawing.")
else
  @assert resb.admissible && length(resb.positions) == n_vertices(Gb.g)
  for e in edges(Gb.g)
    @assert edge_lambda(Gb, resb.positions, e) > 0
  end
  println("Bennis: positive drawing found; convex = $(resb.convex).")
end

println("\n=== All positive-drawing tests passed! ===")
