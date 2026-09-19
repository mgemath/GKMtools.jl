function twisted_matrix_small_quantum_product(V::AbstractGKMVectorBundle, class::GKMClass;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

  G = baseof(V)
  isnothing(beta) || _check_curve_class(G, beta)
  euler_of_V = chern_class(V, rank(V))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the base graph",
  ))
  symmetric_indices = get_symmetric_indices(n)

  pairing = twisted_intersection_matrix(basis, euler_of_V, equivariant = equivariant)

  if !isnothing(beta) && iszero(beta)
    return _MP_inv(pairing) * twisted_intersection_matrix(basis, euler_of_V * class, equivariant = equivariant)
  end
  fast_mode = !equivariant
  
  twist = virtual_zero_section(V)

  if isnothing(beta)
    # coefficient_ring = fast_mode ? QQ : parent(first(localize(class).restrictions))
    # invariant_sums = zero_matrix(coefficient_ring, n, n)
    invariant_sums = twisted_intersection_matrix(basis, euler_of_V * class, equivariant = equivariant) # case iszero(beta)
    for degree in 0:(2 * valency(G))
      for curve in _effective_classes_with_functional_value(GKM_second_homology(G), class, degree)
        iszero(curve) && continue
        show_bar && println("Computing twisted quantum product in curve class $curve:")
        
        invariants = _twisted_invariants(class, basis, symmetric_indices, G, curve, show_bar, equivariant, twist, V)
        
        for ((j, k), invariant) in zip(symmetric_indices, invariants)
          invariant_sums[j, k] += invariant
          j == k || (invariant_sums[k, j] += invariant)
        end
      end
    end

  else
    
    coefficient_ring = fast_mode ? QQ : parent(first(localize(class).restrictions))
    invariant_sums = zero_matrix(coefficient_ring, n, n)

    invariants = _twisted_invariants(class, basis, symmetric_indices, G, beta, show_bar, equivariant, twist, V)
    
    for ((j, k), invariant) in zip(symmetric_indices, invariants)
      invariant_sums[j, k] += invariant
      j == k || (invariant_sums[k, j] += invariant)
    end

  end
  result = _MP_inv(pairing) * invariant_sums

  return result 
end

function twisted_matrix_small_quantum_product_tangent_class(V::AbstractGKMVectorBundle;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)
  
  G = baseof(V)
  @req _test_twisting_positivity_constraint(V) "c1(G)-c1(V) must be strictly positive on every edge"
  class = first_chern_class(G) - first_chern_class(V)
  return twisted_matrix_small_quantum_product(V, class;
    basis = basis, beta = beta, show_bar = show_bar, equivariant = equivariant)
end

function _twisted_invariants(class, basis, symmetric_indices, G, curve, show_bar, equivariant, twist, V)
  fast_mode = !equivariant

  known_zero = fast_mode ?
    _twisted_c1_known_zero(V, curve, class, basis, symmetric_indices) : nothing
  if !isnothing(known_zero) && all(known_zero)
    n = length(basis)
    return zero_matrix(QQ, n, n)
  end

  # class_products = [GKMClass[class, basis[j]] for j in 1:n for k in 1:n]
  marked_insertions = [
    ev(1, class) * ev(2, basis[j]) * ev(3, basis[k]) * twist
    for (j, k) in symmetric_indices
  ]
  invariants = _gromov_witten_gen_0(
    G, curve, 3, marked_insertions, Val(fast_mode);
    show_bar=show_bar, check_degrees=false, g=0, known_zero,
  )
  # marked_insertions = [
  #   ev(1, basis[k]) * twist for j in 1:n for k in 1:n
  # ]
  # invariants = gromov_witten_nomarks(
  #   G, beta, class_products, 1, marked_insertions;
  #   show_bar=show_progress, fast_mode,
  # )
  return invariants 
end

# The Euler twist has degree rank(V) + integral_beta c1(V). Subtracting
# dim Mbar_0,3(G, beta) leaves the degree below, without expanding any
# localization weights. Unknown/nonhomogeneous degrees use the usual check.
function _twisted_c1_known_zero(V, beta, class, basis, indices)
  class_degree = _homogeneous_degree(class)
  basis_degrees = map(_homogeneous_degree, basis)
  (isnothing(class_degree) || any(isnothing, basis_degrees)) && return nothing
  G = baseof(V)
  H2 = GKM_second_homology(G)
  lift = preimage(H2.quotient, beta)
  curve_degree = sum(integrate(class, e) * lift[H2.edge_to_gen[e]] for e in edges(G))
  denominator(curve_degree) == 1 || return nothing
  is_constant(numerator(curve_degree)) || return nothing
  d = Int(constant_coefficient(numerator(curve_degree)))
  return [class_degree + basis_degrees[j] + basis_degrees[k] + rank(V) -
          valency(G) - d != 0 for (j, k) in indices]
end

function _test_twisting_positivity_constraint(V::AbstractGKMVectorBundle)
  class = first_chern_class(baseof(V)) - first_chern_class(V)
  all(e -> begin d = integrate(class, e); denominator(d) == 1 && is_constant(numerator(d)) && constant_coefficient(numerator(d)) > 0 end, edges(baseof(V)))
end
