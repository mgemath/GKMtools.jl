function gkm_graph_of_toric(v::Union{AffineNormalToricVariety, NormalToricVariety}; small_torus::Bool=false)
  @req is_smooth(v) "toric variety must be smooth"
  
  G, labels, flag_rays, edge_flags = _toric_comb_data(v)

  max_cones = maximal_cones(v)
  len = length(max_cones)
  rank = small_torus ? dim(v) : n_rays(v)
  M = free_module(ZZ, rank)
  R = ZZRingElem
  flags::Vector{Vector{ToricFlagWeight{R}}} = Vector{Vector{ToricFlagWeight{R}}}(undef, len)

  for sigma in 1:len

    flags[sigma] = Vector{ToricFlagWeight{R}}(undef, length(rays(max_cones[sigma])))

    for ray_number in 1:length(rays(maximal_cones(v)[sigma]))
      # Store the flag isotropy data for this flag
      # Compute the flag weight using the _omega function
      flags[sigma][ray_number] = ToricFlagWeight{R}(
        _omega(v, sigma, ray_number, M; small_torus=small_torus)
      )

    end
  end
  # return GKMCombinatorialData{ZZRingElem, ToricVertex, ToricFlagWeight{ZZRingElem}}(G, M, labels, flags, edge_flags)
  GKMGraph{ZZRingElem, ToricVertex, ToricFlagWeight{ZZRingElem}}(GKMCombinatorialData{ZZRingElem, ToricVertex, ToricFlagWeight{ZZRingElem}}(G, M, labels, flags, edge_flags), nothing, nothing, nothing, nothing)
end

function gkm_graph_of_orbifold_toric(v::Union{AffineNormalToricVariety, CyclicQuotientSingularity, NormalToricVariety}; small_torus::Bool=false)
  @req is_orbifold(v) "toric variety must be an orbifold"
  G, labels, flag_rays, edge_flags = _toric_comb_data(v)
  M, flags, vertex_isotropy, flag_isotropy = _orbifold_isotropy(v, flag_rays; small_torus=small_torus, R=ZZRingElem)
  return OrbifoldGKMGraph{ZZRingElem, ToricVertex, OrbifoldToricFlagWeight{ZZRingElem}}(GKMCombinatorialData{ZZRingElem, ToricVertex, OrbifoldToricFlagWeight{ZZRingElem}}(G, M, labels, flags, edge_flags), vertex_isotropy, flag_isotropy)
end

function gkm_graph_of_orbifold_toric(v::T; small_torus::Bool=false) where {T <: Union{AbstractStackyFan, AbstractStackyCone}}
  
  G, labels, flag_rays, edge_flags = _toric_comb_data(v, oscar_type=false)
  M, flags, vertex_isotropy, flag_isotropy = _orbifold_isotropy(v, flag_rays; small_torus=small_torus, R=ZZRingElem)
  return OrbifoldGKMGraph{ZZRingElem, ToricVertex, OrbifoldToricFlagWeight{ZZRingElem}}(GKMCombinatorialData{ZZRingElem, ToricVertex, OrbifoldToricFlagWeight{ZZRingElem}}(G, M, labels, flags, edge_flags), vertex_isotropy, flag_isotropy)
end

##### Combinatorial data for stacky fans and cones
function _graph_toric(v::T) where {T <: Union{AffineNormalToricVariety, CyclicQuotientSingularity, NormalToricVariety, AbstractStackyFan, AbstractStackyCone}}
  # Extract the combinatorial data from the toric variety

  max_cones = maximal_cones(v)
  len = length(max_cones)
  g = Graph{Undirected}(len)

  # Build the graph with edges for shared codimension-1 faces
  for sigma1 in 1:(len - 1)
    for sigma2 in (sigma1 + 1):len
      sigma1_cone = max_cones[sigma1]
      sigma2_cone = max_cones[sigma2]

      # Check if these cones share a codimension-1 face
      count(x -> x in rays(sigma1_cone), rays(sigma2_cone)) != (dim(v) - 1) && continue

      add_edge!(g, sigma1, sigma2)
    end
  end
  return g
end

function _toric_comb_data(v; oscar_type::Bool = true)

  num_rays = n_rays(v)
  max_cones = maximal_cones(v)
  len = length(max_cones)
  # number of maximal cones
  g = _graph_toric(v)

  # The labels for the vertices can be the indices of the maximal cones
  labels = ToricVertex.(string.(1:len))

  # The flags and isotropy data will be extracted from the incidence matrix of the cones
  flag_rays::Vector{Vector{Vector{ZZRingElem}}} = Vector{Vector{Vector{ZZRingElem}}}(undef, len)

  for sigma in 1:len
    if oscar_type
      flag_rays[sigma] = Vector{Vector{ZZRingElem}}(undef, length(rays(max_cones[sigma])))
      for (i, r) in enumerate(rays(max_cones[sigma]))
        l = lcm(denominator.(r)...)
        flag_rays[sigma][i] = ZZ.(l .* r)
      end
    else
      flag_rays[sigma] = rays(max_cones[sigma])
    end
  end

  # 4. Build the adjacency dictionary edge_flags between maximal cones
  edge_flags = Dict{Edge,Tuple{Int,Int}}()
  for e in edges(g)
    v1, v2 = src(e), dst(e)
    id_v1 = findfirst(x -> !(x in flag_rays[v2]), flag_rays[v1]) # Find the first common ray between v1 and v2
    id_v2 = findfirst(x -> !(x in flag_rays[v1]), flag_rays[v2]) # Find the first common ray between v1 and v2
    @assert !isnothing(id_v1) "Error: No common ray found between cones $v1 and $v2"
    @assert !isnothing(id_v2) "Error: No common ray found between cones $v1 and $v2"
    edge_flags[e] = (id_v1, id_v2)
  end

  return g, labels, flag_rays, edge_flags
end

function _orbifold_isotropy(v::T, flag_rays::Vector{Vector{Vector{ZZRingElem}}}; small_torus=false, R = ZZRingElem) where {T <: Union{AffineNormalToricVariety, CyclicQuotientSingularity, NormalToricVariety, AbstractStackyFan, AbstractStackyCone}}
  
  max_cones = maximal_cones(v)
  len = length(max_cones)
  rank = small_torus ? dim(v) : n_rays(v)
  M = free_module(ZZ, rank)
  
  flags::Vector{Vector{OrbifoldToricFlagWeight{R}}} = Vector{Vector{OrbifoldToricFlagWeight{R}}}(undef, len)
  vertex_isotropy = Vector{OrbifoldVertexIsotropy}(undef, len)
  flag_isotropy = Vector{Vector{OrbifoldFlagIsotropy}}(undef, len)
  d = dim(v)

  for sigma in 1:len
    invariants, W, U = _invariats_and_weights(flag_rays[sigma], d); println("Invariants for cone $sigma: $invariants", " with weights $W and U $U")

    vertex_isotropy[sigma] = OrbifoldVertexIsotropy(invariants, W)

    flag_isotropy[sigma] = Vector{OrbifoldFlagIsotropy}(
      undef, length(rays(max_cones[sigma]))
    )
    flags[sigma] = Vector{OrbifoldToricFlagWeight{R}}(undef, length(rays(max_cones[sigma])))

    for ray_number in 1:length(rays(maximal_cones(v)[sigma]))
      indices = [i + Int(ray_number <= i) for i in 1:(d - 1)]
      if is_empty(indices)
        # This can happen if the cone is 1-dimensional, in which case we have no rays to include
        flag_isotropy[sigma][ray_number] = smooth_orbifold_flag_isotropy_group(d, 0)
      else

        included_rays = flag_rays[sigma][indices]
        _invariants, _W, _U = _invariats_and_weights(included_rays, d - 1)

        if is_empty(_invariants)
          # This can happen if the flag is smooth, in which case we have no isotropy
          flag_isotropy[sigma][ray_number] = smooth_orbifold_flag_isotropy_group(d, 0)
        else

          # embedding_matrix = U * inv(_U)
          # embedding_matrix = matrix(ZZ, U * inv(_U)[1:2, 1:n_rows(_U)])
          embedding_matrix = sub(U * inv(_U), 1:1, 1:n_rows(_U))

          flag_isotropy[sigma][ray_number] = OrbifoldFlagIsotropy(_invariants, _W, embedding_matrix)
        end
      end
      # Store the flag isotropy data for this flag
      # Compute the flag weight using the _omega function
      flags[sigma][ray_number] = OrbifoldToricFlagWeight{R}(
        _omega(v, sigma, ray_number, M; small_torus=small_torus),
        order_of_generic_stabilizer(
          vertex_isotropy[sigma], flag_isotropy[sigma][ray_number]
        ),
      )

    end
  end
  return M, flags, vertex_isotropy, flag_isotropy
end

function _invariats_and_weights(_rays, d)


  # Edge case: 1-dimensional cone
  M = matrix(ZZ, hcat(_rays...))
  # Step 3: SNF
  S, U, V = snf_with_transform(M)

  # Step 4: invariants
  invariants = [Int(S[i, i]) for i in 1:d if S[i, i] > 1]
  r = length(invariants)

  # Edge case: smooth point
  if r == 0
    return Int[], zeros(Int, 0, Int(d)), identity_matrix(ZZ, Int(n_rows(M)))
  end

  shift = count(i -> S[i, i] == 1, 1:d)

  # # Step 3: correct representation extraction

  W = zeros(Int, r, Int(d))

  for k in 1:r
    for j in 1:d
      # W[k, j] = mod(Int(V[d - j + 1, d - r - shift + k]), invariants[k])
      W[k, j] = mod(Int(V[j, shift + k]), invariants[k])
    end
  end

  return invariants, W, U
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
      den = denominator.(pol_ray)
      lcm_number = length(den) > 1 ? lcm(den...) : den[1]
      ud = lcm_number * pol_ray
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
      _vi = lcm(denominator.(vi)) * vi
      ans += converter(dot(_vi, ud)) * scalars[k]
    end
  end

  return ans
end