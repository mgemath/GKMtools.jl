@doc raw"""
    twisted_c1_matrix(V::AbstractGKMVectorBundle, beta::CurveClass;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the ``q^\beta`` part of the matrix of the `V`-twisted equivariant
quantum product with ``c_1(T_X)-c_1(V)`` on ``X=\operatorname{baseof}(V)``.
The point classes are used by default; pass `basis` to choose another basis.

The twisting class is [`virtual_zero_section`](@ref). The resulting
quantum cohomology maps to that of the smooth zero locus ``Y\subset X`` of a
section of `V`, and ``c_1(T_X)-c_1(V)`` restricts to ``c_1(T_Y)``.

# Arguments
- `V`: a convex GKM vector bundle.
- `beta`: a curve class on `baseof(V)`; the zero class is allowed.
- `basis`: an optional cohomology basis on `baseof(V)`.
- `show_progress`: display progress for the Gromov--Witten computations.
- `fast_mode`: use fast integration and set equivariant parameters to zero.

# Example
```jldoctest twisted_c1_matrix
julia> V = vector_bundle_O(2, [1]);

julia> P2 = baseof(V);

julia> beta0 = zero(GKM_second_homology(P2).H2);

julia> twisted_c1_matrix(V, beta0)
[t1 - t2 - t3                0                0]
[             0   -t1 + t2 - t3                0]
[             0                0   -t1 - t2 + t3]
```
"""
function twisted_c1_matrix(V::AbstractGKMVectorBundle, beta::CurveClass;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  G = baseof(V)
  _check_curve_class(G, beta)
  class = first_chern_class(G) - first_chern_class(V) #tangent class
  euler_of_V = chern_class(V, rank(V))
#   coefficient_ring = parent(first(localize(class).restrictions))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the base graph",
  ))
  symmetric_indices = get_symmetric_indices(n)

  pairing = twisted_intersection_matrix(basis, euler_of_V, equivariant = !fast_mode)

  if iszero(beta)
    return _MP_inv(pairing) * twisted_intersection_matrix(basis, euler_of_V * class, equivariant = !fast_mode)
  end

  known_zero = fast_mode ?
    _twisted_c1_known_zero(V, beta, class, basis, symmetric_indices) : nothing
  if !isnothing(known_zero) && all(known_zero)
    return zero_matrix(QQ, n, n)
  end
#   pairing_2 = matrix(coefficient_ring, n, n, [
#     integrate(euler_of_V * basis[j] * basis[k]) for j in 1:n for k in 1:n
#   ])
#   pairing = fast_mode ? _specialize_matrix_at_origin(G, pairing_2) : pairing_2
  coefficient_ring = fast_mode ? QQ : parent(first(localize(class).restrictions))
  invariant_sums = zero_matrix(coefficient_ring, n, n)
  twist = virtual_zero_section(V)
#   if iszero(beta)
#     for j in 1:n, k in j:n
#       invariant = integrate(euler_of_V * class * basis[j] * basis[k])
#       invariant_sums[j, k] += invariant
#       j == k || (invariant_sums[k, j] += invariant)
#     end
#   else
    # class_products = [GKMClass[class, basis[j]] for j in 1:n for k in 1:n]
    marked_insertions = [
      ev(1, class) * ev(2, basis[j]) * ev(3, basis[k]) * twist
      for (j, k) in symmetric_indices
    ]
    invariants = _gromov_witten_gen_0(
      G, beta, 3, marked_insertions, Val(fast_mode);
      show_bar=show_progress, check_degrees=false, g=0, known_zero,
    )
    # marked_insertions = [
    #   ev(1, basis[k]) * twist for j in 1:n for k in 1:n
    # ]
    # invariants = gromov_witten_nomarks(
    #   G, beta, class_products, 1, marked_insertions;
    #   show_bar=show_progress, fast_mode,
    # )
    for ((j, k), invariant) in zip(symmetric_indices, invariants)
      invariant_sums[j, k] += invariant
      j == k || (invariant_sums[k, j] += invariant)
    end
#   end
  # result = transpose(solve(transpose(pairing), invariant_sums))
  result = _MP_inv(pairing) * invariant_sums
#   fast_mode ? _specialize_matrix_at_origin(G, result) : result
  return result 
end

@doc raw"""
    twisted_c1_matrix_at_q1(V::AbstractGKMVectorBundle;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the matrix of the `V`-twisted equivariant quantum product with
``c_1(T_X)-c_1(V)`` at ``q=1``.

The function sums [`twisted_c1_matrix`](@ref) over all contributing effective
curve classes. It requires ``c_1(T_X)-c_1(V)`` to be strictly positive on every
invariant edge, ensuring that the sum at ``q=1`` is finite.

The point classes are used by default. Pass `basis` to choose another basis.
With `fast_mode=true`, fast integration is used and the returned matrix is
specialized at zero in all equivariant parameters.

Without `fast_mode`, the entries can contain rational functions in the
equivariant parameters. With `fast_mode=true`, the result is a matrix over
`QQ`; the chosen basis must remain a basis after setting those parameters to
zero.

# Example
```jldoctest twisted_c1_matrix_at_q1
julia> V = vector_bundle_O(1, [1]);

julia> M = twisted_c1_matrix_at_q1(V);

julia> size(M)
(2, 2)
```
"""
function twisted_c1_matrix_at_q1(V::AbstractGKMVectorBundle;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  G = baseof(V)
  @req _test_twisting_positivity_constraint(V) "c1(G)-c1(V) must be strictly positive on every edge"
  class = first_chern_class(G) - first_chern_class(V)
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  coefficient_ring = fast_mode ? QQ : parent(first(localize(class).restrictions))
  result = zero_matrix(coefficient_ring, n, n)
  for degree in 0:(2 * valency(G))
    for beta in _effective_classes_with_functional_value(GKM_second_homology(G), class, degree)
      # println(beta)
      # iszero(beta) || continue
      show_progress && println("Computing twisted c1 matrix in curve class $beta:")
      result += twisted_c1_matrix(V, beta; basis, show_progress, fast_mode)
    end
  end
  return result
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
