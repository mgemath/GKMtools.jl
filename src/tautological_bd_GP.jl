export tautological_bd, rank_of_bd

# --- Utility Functions ---

function _check_consistency(R::RootSystem, indices_of_S)
  @req all(i -> 1 <= i <= rank(R), indices_of_S) "Invalid indices"
  return nothing
end

function _same_root_system(lambdas::AbstractArray{WeightLatticeElem})
  isempty(lambdas) && return nothing
  R = root_system(first(lambdas))
  # Check all elements match the root system of the first element
  @req all(l -> root_system(l) == R, lambdas) "Root systems are not the same"
  return nothing
end

# --- Logic ---

function _levi_subroot_system(R::RootSystem, indices_of_S)
  isempty(indices_of_S) && return R
  # Using slicing to preserve Matrix type if necessary, though Oscar usually handles views well.
  cartan_submatrix = cartan_matrix(R)[indices_of_S, indices_of_S]
  levi_root_system = root_system(cartan_submatrix; check=false, detect_type=false)
  return levi_root_system
end

function _levi_subroot_system_and_lambdas_restricted(
  lambdas::AbstractArray{WeightLatticeElem}, indices_of_S
)
  R = root_system(first(lambdas))

  if isempty(indices_of_S)
    # Preserve the shape of lambdas using map
    return R, map(l -> zero(weight_lattice(R)), lambdas)
  end

  levi_root_system = _levi_subroot_system(R, indices_of_S)

  # Pre-calculate range for efficiency
  range_S = 1:length(indices_of_S)

  # Use map to iterate over lambdas, preserving dimensions (Vector, Matrix, etc.)
  lambdas_restricted = map(lambdas) do lambda
    ans = zero(weight_lattice(levi_root_system))
    for index_1 in range_S 
      i = indices_of_S[index_1]
      alpha_i = simple_root(R, i)
      omega_i = fundamental_weight(R, i)
      coe = dot(lambda, alpha_i)//dot(alpha_i, omega_i)
      ans += Int(coe) * fundamental_weight(levi_root_system, index_1)
    end
    ans
    # coeffs = coefficients(lambda)
    # sum(
    #   coeffs[indices_of_S[i]] * dot(lambda, fundamental_weight(R, indices_of_S[i]))//dot(fundamental_weight(levi_root_system, indices_of_S[i]), fundamental_weight(R, indices_of_S[i])) for i in range_S;
    #   init=zero(weight_lattice(levi_root_system)),
    # )
  end

  return levi_root_system, lambdas_restricted
end

# function _lambdas_dynamic(lambdas::AbstractArray{WeightLatticeElem}, indices_of_S)
#   R = root_system(first(lambdas))
#   range_S = 1:length(indices_of_S)

#   # Use map to preserve shape
#   return map(lambdas) do lambda
#     coeffs = coefficients(lambda)
#     # sum(coeffs[indices_of_S[i]] * fundamental_weight(R, indices_of_S[i]) for i in range_S; init=zero(weight_lattice(R)))
#     sum( dot(lambda, ) for i in range_S; init=zero(weight_lattice(R)))
#   end
# end

# function _project(lambda, indices_of_S)
#   R = root_system(lambda)
#   ans = zero_matrix(QQ, rank(R), 1)

#   for i in 1:length(indices_of_S)
#     idx = indices_of_S[i]
#     ans[idx, 1] = dot(lambda, simple_root(R, idx))//dot(simple_root(R, idx), simple_root(R, idx))
#   end

#   return ans
# end

function _lift(R, indices_of_S, lambda_res)
  # This function is not currently used, but it can be helpful for reconstructing full weights from the restricted ones if needed.
  # We would need the original lambda to adjust the lifted weight correctly, so we will pass it as an argument.

  levi_root_system = root_system(lambda_res)
  # coe = coefficients(lambda_res) * transpose(inv(matrix(QQ, cartan_matrix(levi_root_system)))) # need transpose

  emb_coe = zero_matrix(QQ, 1, rank(R))
  range_S = 1:length(indices_of_S)

  for i in range_S
    idx = indices_of_S[i]
    # emb_coe[1, idx] = coe[1, i]
    emb_coe[1, idx] = dot(lambda_res, fundamental_weight(levi_root_system, i))//dot(simple_root(levi_root_system, i), fundamental_weight(levi_root_system, i))
  end
  
  # return emb_coe * cartan_matrix(R)
  return transpose(cartan_matrix(R) * transpose(emb_coe))
end
# function _lift(R, indices_of_S, lambda_res)
#   # This function is not currently used, but it can be helpful for reconstructing full weights from the restricted ones if needed.
#   # We would need the original lambda to adjust the lifted weight correctly, so we will pass it as an argument.

#   levi_root_system = root_system(lambda_res)
#   # coe = coefficients(lambda_res) * transpose(inv(matrix(QQ, cartan_matrix(levi_root_system)))) # need transpose
#   ans = zero(weight_lattice(R))
#   range_S = 1:length(indices_of_S)

#   for i in range_S
#     idx = indices_of_S[i]
#     # emb_coe[1, idx] = coe[1, i]
#     coe = dot(lambda_res, simple_root(levi_root_system, i))//dot(simple_root(levi_root_system, i), fundamental_weight(levi_root_system, i))
#     ans += Int(coe) * fundamental_weight(R, i)
#   end
  
#   # return emb_coe * cartan_matrix(R)
#   return ans
# end

function _rank_of_bd(
  levi_root_system::RootSystem, lambdas_restricted::AbstractArray{WeightLatticeElem}
)
  # map handles multidimensional arrays automatically
  return map(
    l_res -> dim_of_simple_module(Int, levi_root_system, l_res), lambdas_restricted
  )
end

# --- Rank Public Interface ---

function rank_of_bd(lambdas::AbstractArray{WeightLatticeElem}, indices_of_S=Int64[])
  _same_root_system(lambdas)
  R = root_system(first(lambdas))
  _check_consistency(R, indices_of_S)

  levi_root_system, lambdas_restricted = _levi_subroot_system_and_lambdas_restricted(
    lambdas, indices_of_S
  )

  return _rank_of_bd(levi_root_system, lambdas_restricted)
end

function rank_of_bd(lambda::WeightLatticeElem, indices_of_S=Int64[])
  R = root_system(lambda)
  _check_consistency(R, indices_of_S)
  # Wrap in array, call vector version, take first
  return first(rank_of_bd([lambda], indices_of_S))
end

# --- Tautological Bundle Logic ---

function _WP(R, indices_of_S)

  if isempty(indices_of_S)
    return [one(weyl_group(R))] # WP is the whole Weyl group if S is empty
  end

  # No changes to logic, just cleaner iteration
  levi_root_system = _levi_subroot_system(R, indices_of_S)

  return [
    prod(i -> reflection(simple_root(R, indices_of_S[i])), word(a); init=one(weyl_group(R))) for
    a in weyl_group(levi_root_system)
  ]
end

function _cosets_and_repr(Weyl, WP)
  # Optimized lookup using Set, but strictly preserving the Left Coset logic (b * WP)

  total_cosets = div(order(Weyl), length(WP))
  cosets = Vector{Vector{eltype(Weyl)}}(undef, Int(total_cosets))
  reprs = Vector{eltype(Weyl)}(undef, Int(total_cosets))

  # WP is the first coset (Identity * WP)
  cosets[1] = collect(WP)
  reprs[1] = one(Weyl) # Identity

  # Fast lookup to check if element is already covered
  covered = Set{eltype(Weyl)}(WP)

  index = 1

  for b in Weyl
    if b in covered
      continue
    end

    index += 1
    # Original logic: cosets[index] = b .* cosets[index] (where cosets[index] starts as WP)
    # We replicate this: new_coset = b * WP
    new_coset = [b * w for w in WP]

    cosets[index] = new_coset

    # Find representative (shortest length)
    reprs[index] = reduce((x, y) -> length(x) <= length(y) ? x : y, new_coset)

    # Mark these as covered
    union!(covered, new_coset)

    if index == total_cosets
      break
    end
  end

  return cosets, reprs
end

function _gen_matrix_and_type_of_graph(R::RootSystem)
  (fams, ordering) = root_system_type_with_ordering(R)

  type_of_graph = any(fam -> fam[1] in (:E, :F), fams) ? QQFieldElem : ZZRingElem

  return AbstractAlgebra.perm(ordering) *
         block_diagonal_matrix([_generator_matrix(fam) for fam in fams]),
  type_of_graph
end

function _generalized_gkm_flag(R, cosets, reprs, WP, gen_matrix, type_of_graph)
  # Optimized graph construction with O(1) coset lookups

  # Build a lookup for coset index: element -> index
  coset_map = Dict{eltype(reprs),Int}()
  for (idx, coset) in enumerate(cosets)
    for elem in coset
      coset_map[elem] = idx
    end
  end

  WP_set = Set(WP) # Faster lookup

  g = Graph{Undirected}(length(reprs))
  M = free_module(parent(zero(type_of_graph)), n_columns(gen_matrix))
  W = Dict{Edge,AbstractAlgebra.Generic.FreeModuleElem{type_of_graph}}()
  get_root = Dict{Edge,RootSpaceElem}()

  pos_roots = collect(positive_roots(R))
  gens_M = gens(M)
  rank_M = rank(M)

  for i in 1:length(reprs)
    omega = reprs[i]
    inv_omega = inv(omega) # Compute once per outer loop

    for t in pos_roots
      # Remove roots from S
      if reflection(t) in WP_set
        continue
      end

      # Reverted logic: new_rep = omega * reflection(t) (Left multiplication)
      new_rep = omega * reflection(t)

      # Find the representative of the new coset
      # Optimization: Use the dictionary map instead of findfirst loop
      j = get(coset_map, new_rep, nothing)

      if j === nothing
        continue
      end

      # Original Logic: if j > i
      if j > i
        add_edge!(g, j, i) # Orientation: max(i,j) to min(i,j) => j to i

        sign = -1

        # Logic: omega * t = t * inv(omega) (Oscar convention)
        # Reverted logic: vec calculation exactly as original
        vec_coeffs = matrix(parent(zero(type_of_graph)), coefficients(t * inv_omega)*gen_matrix) 

        # Reverted logic: manual sum loop
        val = zero(M)
        for k in 1:rank_M
          val += vec_coeffs[k] * gens_M[k]
        end

        W[Edge(j, i)] = sign * val
      end

      get_root[Edge(i, j)] = t
    end
  end

  ## construct connection
  a = Dict{Tuple{Edge,Edge},ZZRingElem}()

  for _v in vertices(g)
    # Collecting neighbors once is slightly cleaner
    neighbors_v = collect(all_neighbors(g, _v))
    for _w in neighbors_v
      edge_vw = Edge(_v, _w)
      # if !haskey(get_root, edge_vw)
      #   continue
      # end
      alpha = get_root[edge_vw]

      for _u in neighbors_v
        edge_vu = Edge(_v, _u)
        # if !haskey(get_root, edge_vu)
        #   continue
        # end
        beta = get_root[edge_vu]

        a[(edge_vw, edge_vu)] = ZZ(2 * dot(beta, alpha)//dot(alpha, alpha))
      end
    end
  end

  labs = [replace(repr(r), " " => "") for r in reprs]
  GP = gkm_graph(g, labs, M, W)
  con = build_GKM_connection(GP, a)
  set_connection!(GP, con)
  return GP
end

function _tautological_bd(
  lambdas::AbstractArray{WeightLatticeElem}, indices_of_S
)::AbstractArray{GKM_vector_bundle{QQFieldElem}}
  R = root_system(first(lambdas))

  ## Create WP
  Weyl = weyl_group(R)
  WP = _WP(R, indices_of_S)

  # Dynamic Lambdas
  # Using map to preserve shape of 'lambdas'
  levi_root_system, lambdas_restricted = _levi_subroot_system_and_lambdas_restricted(lambdas, indices_of_S)

  # Pre-calculate dominant characters (multiplicities)
  # Optimization: Cache results for identical dynamic lambdas
  dc_m_cache = Dict{WeightLatticeElem,Dict{WeightLatticeElem,Int64}}()

  # We iterate over lambdas to fill the cache, but we don't need a Dict mapping lambda->dc_m anymore
  # if we access it directly in the loop later. However, the original code used a Dict.
  # We will build a direct map that matches the shape of lambdas.
  dc_m_map = map(lambdas_restricted) do ld
    get!(dc_m_cache, ld) do
      dominant_character(levi_root_system, ld)
    end
  end

  cosets, reprs = _cosets_and_repr(Weyl, WP)
  gen_matrix, _ = _gen_matrix_and_type_of_graph(R)

  type_of_graph = QQFieldElem
  # Variety Support
  GP = _generalized_gkm_flag(R, cosets, reprs, WP, gen_matrix, type_of_graph)

  # Pre-calculate Orbits (dc_o)
  # Using map to preserve shape of lambdas
  dc_o_map = Dict{Int64,Dict{WeightLatticeElem,Vector{WeightLatticeElem}}}()

  for i in eachindex(lambdas)
    lambda = lambdas[i]
    lambda_res = lambdas_restricted[i]
    lambda_static = coefficients(lambda) - _lift(R, indices_of_S, lambda_res)#; println("Lambda: ", lambda, " restricted: ", lambda_res, " static part: ", lambda_static, " lifted: ", _lift(R, indices_of_S, lambda_res))
    # lambda_static = lambda - _lift(R, indices_of_S, lambda_res); println("Lambda: ", lambda, " restricted: ", lambda_res, " static part: ", lambda_static, " lifted: ", _lift(R, indices_of_S, lambda_res))
    l_dc_m = dc_m_map[i]
    D = Dict{WeightLatticeElem,Vector{WeightLatticeElem}}()
    for l in keys(l_dc_m)
      l_lifted = _lift(R, indices_of_S, l) # We need to reconstruct the full weight from the restricted one, but we don't have the original lambda here. We can use lambdas[i] and adjust it by the difference between the restricted and original.
      l_int = l_lifted+lambda_static
      dom_weight = sum(k -> Int(l_int[1, k]) * fundamental_weight(R, k), 1:rank(R); init=zero(weight_lattice(R)))
      # dom_weight = l_lifted + lambda_static
      orbit = [dom_weight * w for w in WP]
      D[l] = collect(WeightLatticeElem, Set(orbit)) # Unique elements
    end
    dc_o_map[i] = D
  end
  # dc_o_map = map(lambdas, lambdas_restricted, dc_m_map) do lambda, l_dyn, l_dc_m
  #   D = Dict{WeightLatticeElem,Vector{WeightLatticeElem}}()

  #   for l in keys(l_dc_m)
  #     # Reverted logic: l * w + lambda_static
  #     orbit = [lambda * w for w in WP]; println("Orbit for lambda: ", lambda, " and l: ", l, " is: ", orbit)
  #     D[l] = collect(WeightLatticeElem, Set(orbit)) # Unique elements
  #   end
  #   return D
  # end

  # Sanity Check
  expected_rk = rank_of_bd(lambdas, indices_of_S)

  # Check consistency for each element
  for i in eachindex(lambdas)
    l_dc_o = dc_o_map[i]
    l_dc_m = dc_m_map[i]

    check_val = sum(l -> length(l_dc_o[l]) * l_dc_m[l], keys(l_dc_o); init=0)

    if expected_rk[i] != check_val
      error(
        "Error in tautological vector bundle: computed weights $check_val != expected rank $(expected_rk[i]) at index $i",
      )
    end
  end

  # Construction of Weight Matrices
  M = GP.M
  gens_M = gens(M)
  rank_M = rank(M)
  rank_R = rank(R)

  # Helper function to compute the bundle for a single lambda entry
  # This allows us to use map on the whole array at once
  function compute_single_bundle(index, lambda)
    current_rk = expected_rk[index]
    weightMatrix_S = Matrix{AbstractAlgebra.Generic.FreeModuleElem{type_of_graph}}(
      undef, length(reprs), current_rk
    )

    l_dc_o = dc_o_map[index]
    l_dc_m = dc_m_map[index]

    # Optimization: Sort keys for deterministic column ordering
    # Using height or similar standard sort to ensure stability across runs
    # sorted_keys = sort(collect(keys(l_dc_o)); by=w -> height(w))
    sorted_keys = sort(collect(keys(l_dc_o)); by=w -> sum(coefficients(w)))

    r = 0 # column counter

    for l in sorted_keys #keys(l_dc_o)# 
      orbit_fiber = l_dc_o[l]
      multiplicity = l_dc_m[l]

      for f in orbit_fiber
        # Pre-calculate fiber vector part? No, depends on omega (row)

        # Repeat for multiplicity
        for _ in 1:multiplicity
          r += 1

          # Iterate Rows (Representations)
          for (i, omega) in enumerate(reprs)

            # Reverted logic: t = f * inv(omega)
            t = f * inv(omega)

            # Reverted logic: c * gen_matrix
            # Flattening sum for performance
            # c = coefficients(t) * transpose(inv(matrix(QQ, cartan_matrix(R)))); println("Coefficients for t: ", coefficients(t), " are: ", c)
            # c = coefficients(sum(k -> dot(t, fundamental_weight(R, k)) * simple_root(R, k), 1:rank_R)); println("Coefficients after conversion to simple roots: ", c)
            c = matrix(QQ, 1, rank(R), [dot(t, fundamental_weight(R, k))//dot(simple_root(R, k), fundamental_weight(R, k)) for k in 1:rank_R]); #println("Coefficients after conversion to simple roots, 2nd: ", c)
            vec_coeffs = c * gen_matrix#; println("Vector coefficients for t: ", t, " are: ", vec_coeffs)

            # Reverted logic: -sum(...)
            final_fiber = zero(M)
            for k in 1:rank_M
              final_fiber -= vec_coeffs[k] * gens_M[k]
            end

            weightMatrix_S[i, r] = final_fiber
          end
        end
      end
    end

    M_bd = M
    GMtoM = ModuleHomomorphism(M, M_bd, [gens(M_bd)[i] for i in 1:n_columns(gen_matrix)])

    return vector_bundle(GP, M_bd, GMtoM, weightMatrix_S; calculateConnection=true)
  end

  # Use map to apply computation to every element of lambdas, preserving shape
  # We use CartesianIndices to index into the pre-calculated maps safely
  return map(eachindex(lambdas)) do I
    compute_single_bundle(I, lambdas[I])
  end
end

function tautological_bd(
  lambda::WeightLatticeElem, indices_of_S=Int64[]
)::GKM_vector_bundle{QQFieldElem}
  R = root_system(lambda)
  _check_consistency(R, indices_of_S)
  return first(tautological_bd([lambda], indices_of_S))
end

function tautological_bd(
  lambdas::AbstractArray{WeightLatticeElem}, indices_of_S=Int64[]
)::AbstractArray{GKM_vector_bundle{QQFieldElem}}
  _same_root_system(lambdas)
  R = root_system(first(lambdas))
  _check_consistency(R, indices_of_S)

  return _tautological_bd(lambdas, indices_of_S)
end