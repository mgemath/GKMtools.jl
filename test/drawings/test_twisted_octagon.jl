using GKMtools
using Oscar

# Twisted 8-gon with weights x, y, -x, -y, x, y, -x, -y
println("=== Twisted 8-gon ===")
w = [1 0; 0 1; -1 0; 0 -1; 1 0; 0 1; -1 0; 0 -1]
G = gkm_2d(w)
println("Vertices: $(n_vertices(G.g)), rank: $(rank_torus(G))")

basis, l_e = drawing_space(G)
r = ncols(basis)
println("dim D(G) = $r")

println("\n--- Without vertex injectivity ---")
reps = admissible_drawing_representatives(G)
println("Admissible components: $(length(reps))")

# Check which drawings have vertex collisions
for (i, pos) in enumerate(reps)
  n = length(pos)
  collisions = Tuple{Int,Int}[]
  for v in 1:n, w in (v+1):n
    if pos[v] == pos[w]
      push!(collisions, (v, w))
    end
  end
  if !isempty(collisions)
    println("  Rep $i has collisions: $collisions")
  else
    println("  Rep $i is vertex-injective")
  end
end

println("\n--- With vertex injectivity ---")
reps_vi = admissible_drawing_representatives(G; require_vertex_injectivity=true)
println("Vertex-admissible components: $(length(reps_vi))")

for (i, pos) in enumerate(reps_vi)
  n = length(pos)
  for v in 1:n, w in (v+1):n
    @assert pos[v] != pos[w] "Rep $i: vertices $v and $w collide!"
  end
  println("  Rep $i is vertex-injective ✓")
end

println("\n=== Done ===")
