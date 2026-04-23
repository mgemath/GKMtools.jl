using Oscar, GKMtools

# ─── helpers ──────────────────────────────────────────────────────────────────

function verify_drawing(G, positions)
  d = rank_torus(G)
  for e in edges(G.g)
    v, w = src(e), dst(e)
    we = [QQ(GKMtools._w(G, e)[k]) for k in 1:d]
    delta = positions[w] - positions[v]
    # Check that delta is proportional to we: delta[j]*we[k] == delta[k]*we[j] for all j < k
    for j in 1:d, k in j+1:d
      @assert delta[j] * we[k] == delta[k] * we[j] "Drawing invalid for edge $e at (j=$j, k=$k)"
    end
    # Check admissibility: delta != 0
    @assert any(c -> c != 0, delta) "Drawing not admissible: edge $e collapses"
  end
  return true
end

# ─── Test 1: P^2 (flag_variety([2,1])) ────────────────────────────────────────
println("=== Test 1: P^2 = flag_variety([2,1]) ===")
G = flag_variety(GKM_graph, [2, 1])
println("n_vertices = $(n_vertices(G.g)), rank_torus = $(rank_torus(G)), valency = $(valency(G))")

basis, l_e = drawing_space(G)
println("dim D(G) = $(ncols(basis))")
@assert ncols(basis) == 1 "Expected 1-dimensional drawing space for P^2"

reps = admissible_drawing_representatives(G)
println("Number of admissible components (up to sign) = $(length(reps))")
@assert length(reps) == 1 "Expected 1 admissible component for P^2"

# Verify the representative is a valid admissible drawing
verify_drawing(G, reps[1])
println("Representative verified as valid admissible drawing.")
println("Vertex positions:")
for (i, p) in enumerate(reps[1])
  println("  vertex $i => $p")
end

# ─── Test 2: P^3 (flag_variety([1,3])) ────────────────────────────────────────
println("\n=== Test 2: P^3 = flag_variety([1,3]) ===")
G3 = flag_variety(GKM_graph, [1, 3])
println("n_vertices = $(n_vertices(G3.g)), rank_torus = $(rank_torus(G3)), valency = $(valency(G3))")

basis3, l_e3 = drawing_space(G3)
println("dim D(G) = $(ncols(basis3))")

reps3 = admissible_drawing_representatives(G3)
println("Number of admissible components (up to sign) = $(length(reps3))")
@assert length(reps3) >= 1 "Expected at least 1 admissible component"

for (i, pos) in enumerate(reps3)
  verify_drawing(G3, pos)
end
println("All representatives verified as valid admissible drawings.")

# ─── Test 3: Hirzebruch surface (flag_variety([1,1,1,1]) is the full flag) ───
println("\n=== Test 3: Full flag variety F(1,1,1) = flag_variety([1,1,1]) ===")
F3 = flag_variety(GKM_graph, [1, 1, 1])
println("n_vertices = $(n_vertices(F3.g)), rank_torus = $(rank_torus(F3)), valency = $(valency(F3))")

basisF3, l_eF3 = drawing_space(F3)
println("dim D(G) = $(ncols(basisF3))")

repsF3 = admissible_drawing_representatives(F3)
println("Number of admissible components (up to sign) = $(length(repsF3))")
@assert length(repsF3) >= 1 "Expected at least 1 admissible component"

for pos in repsF3
  verify_drawing(F3, pos)
end
println("All representatives verified.")

# ─── Test 4: projective_space ─────────────────────────────────────────────────
println("\n=== Test 4: projective_space(GKM_graph, 2) ===")
P2 = projective_space(GKM_graph, 2)
println("n_vertices = $(n_vertices(P2.g)), rank_torus = $(rank_torus(P2)), valency = $(valency(P2))")

basisP2, _ = drawing_space(P2)
println("dim D(G) = $(ncols(basisP2))")

repsP2 = admissible_drawing_representatives(P2)
println("Number of admissible components (up to sign) = $(length(repsP2))")
for pos in repsP2
  verify_drawing(P2, pos)
end
println("All representatives verified.")


# ─── Test 5: project_drawings ─────────────────────────────────────────────────
println("\n=== Test 5: project_drawings ===")
G = flag_variety(GKM_graph, [2, 1])   # d = 3
reps = admissible_drawing_representatives(G)

# Default (random) projection
proj = project_drawings(reps)
println("Default projection: $(length(proj)) drawing(s), each has $(length(proj[1])) vertices")
@assert length(proj) == length(reps)
@assert all(p -> p isa Tuple{Float64,Float64}, proj[1])

# Explicit 2×3 projection
P_explicit = zero_matrix(QQ, 2, 3)
P_explicit[1, 1] = QQ(1); P_explicit[1, 3] = QQ(-1)
P_explicit[2, 2] = QQ(1); P_explicit[2, 3] = QQ(-1)
proj2 = project_drawings(reps; projection = P_explicit)
println("Explicit projection, vertex 2 position: $(proj2[1][2])")
@assert proj2[1][1] == (0.0, 0.0)   # vertex 1 always at origin

# flag_variety([1,2]) has d=3; use projective_space which also has d=3
# Test with a d=2 graph: flag_variety([1,1]) is P^1, d=2
G_2d = flag_variety(GKM_graph, [1, 1])
reps_2d = admissible_drawing_representatives(G_2d)
println("P^1 (d=2): $(length(reps_2d)) component(s)")
proj_2d = project_drawings(reps_2d)
println("Identity projection on P^1, vertex positions: $(proj_2d[1])")


# ─── Test 6: invariance under GKM isomorphism (regression for find_initial_sv) ─
# The first almost-positive survivor is GKM-isomorphic to gkm_3d_twisted_flag().
# Previously, the two graphs gave different numbers of admissible drawing
# representatives (5 vs 0) because the initial sign-vector search only tried
# ±e_k and ±1 points, all of which happened to lie on a hyperplane of the
# arrangement in the survivor's basis.  The moment-curve probe now used in
# find_initial_sv() is basis-independent, so counts must agree.
println("\n=== Test 6: twisted flag vs its isomorphic survivor ===")
F = gkm_3d_twisted_flag()

# Build the survivor G inline from the edge weights listed in the bug report
# (avoids any dependency on the ap_classification/ data pipeline).
g6 = Graph{Undirected}(6)
add_edge!(g6, 2, 1)
add_edge!(g6, 3, 1)
add_edge!(g6, 4, 1)
add_edge!(g6, 5, 2)
add_edge!(g6, 5, 3)
add_edge!(g6, 5, 4)
add_edge!(g6, 6, 2)
add_edge!(g6, 6, 3)
add_edge!(g6, 6, 4)
M6 = free_module(ZZ, 2)
e1, e2 = gens(M6)
w6 = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(
  Edge(2, 1) =>   e1 + e2,
  Edge(3, 1) =>   e1,
  Edge(4, 1) =>      - e2,
  Edge(5, 2) =>   e1,
  Edge(5, 3) =>      - e2,
  Edge(5, 4) =>   e1 - e2,
  Edge(6, 2) =>      - e2,
  Edge(6, 3) =>   e1 - e2,
  Edge(6, 4) =>   e1 - 2*e2,
)
G = gkm_graph(g6, ["v1","v2","v3","v4","v5","v6"], M6, w6)

reps_F = admissible_drawing_representatives(F)
reps_G = admissible_drawing_representatives(G)
println("F gave $(length(reps_F)) components, G gave $(length(reps_G)).")
@assert length(reps_F) == length(reps_G) "admissible component count must be GKM-iso invariant"
@assert length(reps_F) == 5 "expected 5 admissible components for the twisted-flag/survivor graph"
for pos in reps_G
  verify_drawing(G, pos)
end
println("Invariance test passed: both F and G give $(length(reps_F)) components.")

println("\n=== All tests passed! ===")
