using GKMtools
using Oscar

# Test 1: Basic cases — flag varieties should have the same number of components
# with and without vertex injectivity (all vertices are already distinct)
println("=== Test 1: flag_variety([2,1]) ===")
G1 = flag_variety(GKM_graph, [2, 1])
reps_normal = admissible_drawing_representatives(G1)
reps_vi = admissible_drawing_representatives(G1; require_vertex_injectivity=true)
println("Without vertex injectivity: $(length(reps_normal)) components")
println("With vertex injectivity: $(length(reps_vi)) components")
@assert length(reps_normal) == length(reps_vi) "CP2: counts should match"

println("\n=== Test 2: flag_variety([1,1,1]) ===")
G2 = flag_variety(GKM_graph, [1, 1, 1])
reps_normal2 = admissible_drawing_representatives(G2)
reps_vi2 = admissible_drawing_representatives(G2; require_vertex_injectivity=true)
println("Without vertex injectivity: $(length(reps_normal2)) components")
println("With vertex injectivity: $(length(reps_vi2)) components")
@assert length(reps_normal2) == length(reps_vi2) "F(1,1,1): counts should match"

# Check that all returned drawings are actually vertex-injective
for (i, pos) in enumerate(reps_vi2)
  n = length(pos)
  for v in 1:n, w in (v+1):n
    @assert pos[v] != pos[w] "F(1,1,1) rep $i: vertices $v and $w collide!"
  end
end
println("All F(1,1,1) drawings are vertex-injective ✓")

# Test 3: The custom example from vertex_injectivity.jl where vertices may collide
println("\n=== Test 3: Custom graph with potential vertex collisions ===")
include("vertex_injectivity.jl")
G3 = get_vertex_inj_example()
reps_normal3 = admissible_drawing_representatives(G3)
println("Without vertex injectivity: $(length(reps_normal3)) components")
reps_vi3 = admissible_drawing_representatives(G3; require_vertex_injectivity=true)
println("With vertex injectivity: $(length(reps_vi3)) components")

# Check vertex injectivity of all returned drawings
for (i, pos) in enumerate(reps_vi3)
  n = length(pos)
  for v in 1:n, w in (v+1):n
    @assert pos[v] != pos[w] "Custom graph rep $i: vertices $v and $w collide!"
  end
end
println("All custom graph drawings are vertex-injective ✓")

# Test 4: convex_drawing_representatives passes through the kwarg
println("\n=== Test 4: convex_drawing_representatives with vertex injectivity ===")
creps_vi = convex_drawing_representatives(G2; require_vertex_injectivity=true)
println("F(1,1,1) convex + vertex-injective: $(length(creps_vi)) components")

println("\n=== All tests passed! ===")
