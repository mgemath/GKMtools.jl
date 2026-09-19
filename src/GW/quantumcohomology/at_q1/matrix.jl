function matrix_small_quantum_product(G::AbstractGKMGraph, class::GKMClass;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

  @req all(e -> chern_number(e, G) > 0, edges(G)) "the graph must be strictly nef"
  class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the given graph",
  ))

  pairing = intersection_matrix(basis; equivariant = equivariant)

  if !isnothing(beta) && iszero(beta)
    return _MP_inv(pairing) * twisted_intersection_matrix(basis, class, equivariant = equivariant)
  end

  symmetric_indices = get_symmetric_indices(n)

  if isnothing(beta)
    invariant_sums = twisted_intersection_matrix(basis, class; equivariant = equivariant) # case iszero(curve)
    for degree in 0:(2 * valency(G))
      for curve in _effective_classes_with_functional_value(
          GKM_second_homology(G), first_chern_class(G), degree)
        if iszero(curve)
          continue
        end
        show_bar && println("Computing quantum product in curve class $curve:")        
        invariants = _invariants(class, basis, symmetric_indices, G, curve, show_bar, equivariant)
        
        for ((j, k), invariant) in zip(symmetric_indices, invariants)
          invariant_sums[j, k] += invariant
          j == k || (invariant_sums[k, j] += invariant)
        end
      end
    end
  else
    invariants = _invariants(class, basis, symmetric_indices, G, beta, show_bar, equivariant)
    coefficient_ring = equivariant ? parent(first(localize(class).restrictions)) : QQ
    invariant_sums = zero_matrix(coefficient_ring, n, n)
        
    for ((j, k), invariant) in zip(symmetric_indices, invariants)
      invariant_sums[j, k] += invariant
      j == k || (invariant_sums[k, j] += invariant)
    end
  end

  result = _MP_inv(pairing) * invariant_sums

  return result
end

function matrix_small_quantum_product_tangent_class(G::AbstractGKMGraph;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)
  
  class = first_chern_class(G)
  return matrix_small_quantum_product(G, class;
    basis = basis, beta = beta, show_bar = show_bar, equivariant = equivariant)
end

function _invariants(class, basis, symmetric_indices, G, curve, show_bar, equivariant)

  fast_mode = !equivariant
  marked_insertions = [
  ev(1, class) * ev(2, basis[j]) * ev(3, basis[k])
  for (j, k) in symmetric_indices
  ]
  invariants = gromov_witten(
    G, curve, 3, marked_insertions; show_bar = show_bar, fast_mode,
  )
  # marked_insertions = [
  #   GKMClass[class, basis[j], basis[k]]
  #   for (j, k) in symmetric_indices
  # ]
  # invariants = gromov_witten_nomarks(
  #   G, curve, marked_insertions; show_bar=show_progress, fast_mode,
  # )
  return invariants
end
