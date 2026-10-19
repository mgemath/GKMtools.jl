@doc raw"""
    twisted_matrix_small_quantum_product(V, class; basis=nothing, beta=nothing, show_bar=true, equivariant=false)
    twisted_matrix_small_quantum_product(V, classes::Vector; kwargs...)

Return the ``q^\beta`` part of the matrix of the `V`-twisted equivariant
quantum product with ``c_1(T_X)-c_1(V)`` on ``X=\operatorname{baseof}(V)``.
The point classes are used by default; pass `basis` to choose another basis.

The twisting class is [`virtual_zero_section`](@ref). The resulting
quantum cohomology maps to that of the smooth zero locus ``Y\subset X`` of a
section of `V`, and ``c_1(T_X)-c_1(V)`` restricts to ``c_1(T_Y)``.

Twisted version of [`matrix_small_quantum_product`](@ref) for GKM vector bundles.

# Example
```jldoctest twisted_c1_matrix
julia> V = vector_bundle_O(5, [3]);

julia> P5 = baseof(V);

julia> H = first_chern_class(V)/3
GKM class with restrictions: 
[1//3*t1, -2//3*t1 + t2, -2//3*t1 + t3, -2//3*t1 + t4, -2//3*t1 + t5, -2//3*t1 + t6]

julia> basis = [H^i for i in 0:(num_vertices(P5)-2)];

julia> M = twisted_matrix_small_quantum_product_tangent_class(V, basis = basis, show_bar = false)
[0   0   18    0    0]
[3   0    0   45    0]
[0   3    0    0   18]
[0   0    3    0    0]
[0   0    0    3    0]
```
"""
function twisted_matrix_small_quantum_product(V::AbstractGKMVectorBundle, class::GKMClass;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

  return only(twisted_matrix_small_quantum_product(V, [class];
    basis, beta, show_bar, equivariant))
end

function twisted_matrix_small_quantum_product(V::AbstractGKMVectorBundle,
    classes::Vector{GKMClass{GT,R}};
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false) where {GT,R}

  G = baseof(V)
  all(class -> class.graph === G, classes) || throw(ArgumentError(
    "all classes must belong to the base graph",
  ))
  isnothing(beta) || _check_curve_class(G, beta)
  isempty(classes) && return []
  euler_of_V = chern_class(V, rank(V))
  # euler_of_V = chern_class(V, 0)
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph == G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the base graph",
  ))
  symmetric_indices = get_symmetric_indices(n)

  inverse_pairing = _MP_inv(twisted_intersection_matrix(basis, euler_of_V; equivariant))

  if !isnothing(beta) && iszero(beta)
    return [inverse_pairing * twisted_intersection_matrix(basis, euler_of_V * class; equivariant)
            for class in classes]
  end
  twist = virtual_zero_section(V)
  # twist = reduced_virtual_zero_section(V)

  if isnothing(beta)
    invariant_sums = [twisted_intersection_matrix(basis, euler_of_V * class; equivariant)
                      for class in classes]
    twisted_c1 = first_chern_class(G) - first_chern_class(V)
    # twisted_c1 = chern_class(G, 1)
    for degree in 0:(2 * valency(G))
      for curve in _effective_classes_with_functional_value(GKM_second_homology(G), twisted_c1, degree)
        iszero(curve) && continue
        show_bar && println("Computing twisted quantum product in curve class $curve:")

        invariants = _twisted_invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant, twist, V)

        for ((i, (j, k)), invariant) in zip(
            ((i, pair) for i in eachindex(classes) for pair in symmetric_indices), invariants)
          invariant_sums[i][j, k] += invariant
          j == k || (invariant_sums[i][k, j] += invariant)
        end
      end
    end

  else

    Z = zero_matrix(
      equivariant ? parent(first(localize(first(classes)).restrictions)) : QQ, n, n,
    )
    invariant_sums = [deepcopy(Z) for class in classes]

    invariants = _twisted_invariants(classes, basis, symmetric_indices, G, beta, show_bar, equivariant, twist, V)

    for ((i, (j, k)), invariant) in zip(
            ((i, pair) for i in eachindex(classes) for pair in symmetric_indices), invariants)
      invariant_sums[i][j, k] += invariant
      j == k || (invariant_sums[i][k, j] += invariant)
    end

  end
  return [inverse_pairing * M for M in invariant_sums]
  # return [ M for M in invariant_sums]
end

function twisted_matrix_small_quantum_product_tangent_class(V::AbstractGKMVectorBundle;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

  G = baseof(V)
  @req _test_twisting_positivity_constraint(V) "c1(G)-c1(V) must be strictly positive on every edge"
  class = first_chern_class(G) - first_chern_class(V)
  return twisted_matrix_small_quantum_product(V, class;
    basis = basis, beta = beta, show_bar = show_bar, equivariant = equivariant)
end

function _twisted_invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant, twist, V)
  fast_mode = !equivariant
  zero_checks = fast_mode ?
    [_twisted_c1_known_zero(V, curve, class, basis, symmetric_indices) for class in classes] : nothing
  known_zero = isnothing(zero_checks) || any(isnothing, zero_checks) ?
    nothing : reduce(vcat, zero_checks)
  if !isnothing(known_zero) && all(known_zero)
    return [QQ(0) for _ in known_zero]
  end

  class_products = [
    GKMClass[class, basis[j], basis[k]]
    for class in classes for (j, k) in symmetric_indices
  ]
  return gromov_witten_nomarks(
    G, curve, class_products, 0, twist; show_bar, fast_mode, known_zero,
  )
  # P = [ev(1, class) * ev(2, basis[j]) * ev(3, basis[k]) * twist
  #      for class in classes for (j, k) in symmetric_indices]
  # return gromov_witten(
  #   G, curve, 3, P; fast_mode = true, show_bar = show_bar
  # )
end

# function _twisted_invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant, twist, V)
#   fast_mode = !equivariant

#   zero_checks = fast_mode ?
#     [_twisted_c1_known_zero(V, curve, class, basis, symmetric_indices) for class in classes] : nothing
#   known_zero = isnothing(zero_checks) || any(isnothing, zero_checks) ?
#     nothing : reduce(vcat, zero_checks)
#   if !isnothing(known_zero) && all(known_zero)
#     return [QQ(0) for _ in known_zero]
#   end

#   marked_insertions = [
#     ev(1, class) * ev(2, basis[j]) * ev(3, basis[k]) * twist
#     for class in classes for (j, k) in symmetric_indices
#   ]
#   invariants = _gromov_witten_gen_0(
#     G, curve, 3, marked_insertions, Val(fast_mode);
#     show_bar=show_bar, check_degrees=false, g=0, known_zero,
#   )
#   return invariants
# end

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
  twisted_c1 = first_chern_class(G) - first_chern_class(V)
  curve_degree = sum(integrate(twisted_c1, e) * lift[H2.edge_to_gen[e]] for e in edges(G))
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
