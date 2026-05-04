@doc raw"""
    gkm_graph_of_toric(v::NormalToricVariety; small_torus::Bool=false) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the smooth toric variety `v`.

If the variety is projective, all flags will be connected to edges.
If the variety is non-projective, codimension-1 faces that are not shared by two maximal cones
will correspond to standalone flags.

# Dimension of the torus
If the optional argument `small_torus` is `false` (default value) then the torus rank of the
result is the number of rays of `v`.
If `small_torus` is `true` then the torus rank of the result is the dimension of `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> gkm_graph_of_toric(P2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1)
3 -> 1 => (0, 1, -1)
3 -> 2 => (-1, 1, 0)

julia> gkm_graph_of_toric(P2; small_torus=true)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0)
3 -> 1 => (0, 1)
3 -> 2 => (-1, 1)

julia> F = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> gkm_graph_of_toric(F)
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (3, 1, 0, -1)
4 -> 1 => (0, 1, 3, -1)
4 -> 3 => (-1, 0, 1, 0)

julia> gkm_graph_of_toric(F; small_torus=true)
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0)
3 -> 2 => (3, 1)
4 -> 1 => (0, 1)
4 -> 3 => (-1, 0)

julia> gkm_graph_of_toric(affine_space(NormalToricVariety, 3))
GKM graph with 1 nodes, valency 3 and axial function:
Standalone flags:
1.1 => (-1, 0, 0)
1.2 => (0, -1, 0)
1.3 => (0, 0, -1)
```
"""

function gkm_graph_of_toric(
  v::NormalToricVariety; small_torus::Bool=false
)::GKMGraph{ZZRingElem}
  return _gkm_graph_of_toric(v; small_torus, check_smothness=true, return_type=ZZRingElem)
end

function _gkm_graph_of_toric(
  v; small_torus::Bool=false, check_smothness::Bool=true,
  return_type::Type{R}
) where {R}
  @req !check_smothness || is_smooth(v) "toric variety must be smooth"

  len = length(maximal_cones(v))
  g = Graph{Undirected}(len)
  base_ring_of_M = parent(zero(R))
  M = free_module(base_ring_of_M, small_torus ? dim(v) : n_rays(v))
  W = Dict{Edge,AbstractAlgebra.Generic.FreeModuleElem{return_type}}()

  # Build the graph with edges for shared codimension-1 faces
  for sigma1 in 1:(len - 1)
    for sigma2 in (sigma1 + 1):len
      sigma1_cone = maximal_cones(v)[sigma1]
      sigma2_cone = maximal_cones(v)[sigma2]

      # Check if these cones share a codimension-1 face
      count(x -> x in rays(sigma1_cone), rays(sigma2_cone)) != (dim(v) - 1) && continue

      # Find the ray in sigma1 that's not in sigma2
      ray1 = findfirst(r -> !(r in rays(sigma2_cone)), rays(sigma1_cone))

      add_edge!(g, sigma1, sigma2)
      W[Edge(sigma2, sigma1)] = _omega(v, sigma1, ray1, M; small_torus)
    end
  end

  # Create the base GKM graph
  G = gkm_graph(g, ["$i" for i in 1:len], M, W; check=false) # don't check as valency is wrong until flags are added.

  # Now add standalone flags for codimension-1 faces that don't belong to another maximal cone
  for sigma_idx in 1:len
    sigma = maximal_cones(v)[sigma_idx]

    # For each ray in this cone
    for (ray_idx, ray) in enumerate(rays(sigma))
      # Check if this ray corresponds to a codimension-1 face shared with another cone
      # by checking if removing this ray gives a face contained in another maximal cone
      codim1_face = [r for r in rays(sigma) if r != ray]

      # Check if this codim-1 face is contained in any other maximal cone
      is_shared = false
      for other_idx in 1:len
        other_idx == sigma_idx && continue

        other_sigma = maximal_cones(v)[other_idx]
        if all(r -> r in rays(other_sigma), codim1_face)
          is_shared = true
          break
        end
      end

      # If not shared, this is a standalone flag
      if !is_shared
        weight = _omega(v, sigma_idx, ray_idx, M; small_torus)
        # add_standalone_flag!(G, sigma_idx, -weight)
        G.flags[sigma_idx] = vcat(G.flags[sigma_idx], GKMFlag(sigma_idx, -weight, nothing))
        # push!(flags[sigma_idx], GKMFlag(sigma_idx, -weight, nothing))
      end
    end
  end

  #   @req isvalid(G) "gkm_graph_of_toric produced invalid result."
  return G
  # return GKMGraph(g, M, ["$i" for i in 1:len], flags, Dict{Edge,Tuple{Int,Int}}(), nothing, nothing, nothing, nothing)
end

###############################################################################
function _omega(v,
  n_SIGMA::Int64,
  ray_idx::Int64,
  M::AbstractAlgebra.Generic.FreeModule{R}; small_torus::Bool=false,
)::AbstractAlgebra.Generic.FreeModuleElem{R} where {R}
  SIGMA = maximal_cones(v)[n_SIGMA]
  ray = rays(SIGMA)[ray_idx]
  scalars = gens(M)
  ud = QQFieldElem[]

  # Find a polarizing ray that pairs non-trivially with our ray
  for pol_ray in rays(polarize(SIGMA))
    if dot(ray, pol_ray) != 0
      ud = lcm(denominator.(pol_ray)) * pol_ray
      break
    end
  end

  ans = zero(M)

  converter(x::QQFieldElem)::R = base_ring(M) === ZZ ? numerator(x) : x

  if small_torus
    for k in 1:dim(v)
      ans += converter(ud[k]) * scalars[k]
    end
  else
    for (k, vi) in enumerate(rays(v))
      ans += converter(dot(vi, ud)) * scalars[k]
    end
  end

  return ans
end

function orbifold_gkm_graph_of_toric(X; small_torus::Bool=false)
  @req is_orbifold(X) "toric variety must be an orbifold"

  base = _gkm_graph_of_toric(X; small_torus, check_smothness=false, return_type=QQFieldElem)

  vertex_isotropy = Vector{VertexIsotropyData}(undef, num_vertices(base))
  flag_isotropy_data = Dict{Flag, FlagIsotropyData}()

  for v in eachindex(vertex_isotropy)
    SIGMA = maximal_cones(X)[v]

    vertex_isotropy[v] = _vertex_isotropy_data(SIGMA)
    # mat = []
    # # mat = matrix(QQ, hcat(rays(SIGMA)...))
    # for r in rays(SIGMA)
    #   rd = lcm(denominator.(r))*r
    #   push!(mat, rd)
    # end
    # mat = hcat(mat...)
    # smith = snf(matrix(ZZ, mat))
    # vertex_isotropy[v] = IsotropyData(Int(prod(smith[i, i] for i in 1:dim(X))), [Int(smith[i, i]) for i in 1:dim(X)])
  end

  for (i, sub_c) in enumerate(cones(X, dim(X) - 1))
    invariants, W = _invariats_and_weights(rays(sub_c), dim(X) - 1)
    f = Flag(i, nothing)
    flag_isotropy_data[f] = FlagIsotropyData(invariants, W)

  end

  # for _v in 1:num_vertices(base)
  #   SIGMA1 = maximal_cones(X)[_v]
  #   for _w in 1:num_vertices(base)
    
  #   SIGMA2 = maximal_cones(X)[_w]
  #   f = Flag(src(e), Edge(_))
  #   println(_flag_isotropy_data(SIGMA1, SIGMA2))
  #   flag_isotropy_data[f] = _flag_isotropy_data(SIGMA1, SIGMA2)
  # end
  # return flag_isotropy_data

  return OrbifoldGKMGraph(base, vertex_isotropy, flag_isotropy_data)
end

function _flag_isotropy_data(SIGMA1, SIGMA2)
  d = length(rays(SIGMA1)) - 1
  face = [r for r in rays(SIGMA1) if r in rays(SIGMA2)]
  invariants, W = _invariats_and_weights(face, d)

  return FlagIsotropyData(invariants, W)
end

function _vertex_isotropy_data(SIGMA)
  d = length(rays(SIGMA))

  invariants, W = _invariats_and_weights(rays(SIGMA), d)

  return VertexIsotropyData(invariants, W)
end

function _invariats_and_weights(_rays, d)
  # Step 1: clear denominators
  cols = Vector{Vector{Int}}(undef, d)
  for (i, r) in enumerate(_rays)
    l = lcm(denominator.(r)...)
    cols[i] = Int.(l .* r)
  end
  # Step 2: matrix
  M = matrix(ZZ, hcat(cols...))


  # Step 3: SNF
  S, U, V = snf_with_transform(M)
  # Step 4: invariants
  invariants = [Int(S[i, i]) for i in 1:d if S[i, i] > 1]
  shift = count(i -> S[i, i] == 0, 1:d)
  if shift > 0
    println("M is not full rank, this behaviour is not tested")
  end
  r = length(invariants)

  # Edge case: smooth point
  if r == 0
    return Int[], zeros(Int, 0, d)
  end

  # # Step 3: correct representation extraction

  W = zeros(Int, r, d)

  for k in 1:r
    for j in 1:d
      W[k, j] = mod(Int(V[d - j + 1, d - r - shift + k]), invariants[k])
    end
  end

  return invariants, W
end
