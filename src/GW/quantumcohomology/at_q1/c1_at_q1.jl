function _effective_classes_with_functional_value(H2::GKM_H2, class::GKMClass,
    value::Integer; check::Bool=true)
  G = class.graph
  integers = free_module(ZZ, 1)
  z = first(gens(integers))
  edge_values = Vector{typeof(z)}(undef, rank(H2.edge_lattice))
  for (e, edge_index) in H2.edge_to_gen
    degree = integrate(class, e)
    @req denominator(degree) == 1 "the class does not have integral degree on $e"
    @req is_constant(numerator(degree)) "the degree on $e is not constant"
    degree = ZZ(constant_coefficient(numerator(degree)))
    check && @req degree > 0 "the class must be strictly positive on every edge"
    edge_values[edge_index] = degree * z
  end
  edge_to_value = ModuleHomomorphism(H2.edge_lattice, integers, edge_values)
  H2_to_value = ModuleHomomorphism(H2.H2, integers,
    [edge_to_value(preimage(H2.quotient, g)) for g in gens(H2.H2)])
  success, beta0 = has_preimage_with_preimage(H2_to_value, value * z)
  success || return CurveClass[]
  K, inclusion = kernel(H2_to_value)
  inclusion_matrix = transpose(matrix(inclusion))
  dual_rays = rays(H2.dual_cone)
  inequalities = matrix(QQ, length(dual_rays), rank(H2.H2),
    [dual_rays[i][j] for i in eachindex(dual_rays) for j in 1:rank(H2.H2)])
  rhs = inequalities * [beta0[i] for i in 1:rank(H2.H2)]
  polytope = polyhedron(-inequalities * inclusion_matrix, rhs)
  (beta0 + inclusion(K([p[i] for i in 1:rank(K)])) for p in lattice_points(polytope))
end

function _multiplication_matrix_at_q1(G::AbstractGKMGraph, class::GKMClass;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  @req all(e -> chern_number(e, G) > 0, edges(G)) "the graph must be strictly nef"
  class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  coefficient_ring = parent(first(localize(class).restrictions))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the given graph",
  ))

  # If P[j, k] = integral(T_j T_k), first collect
  # S[j, k] = sum_beta GW(class, T_j, T_k). The multiplication matrix is
  # then the solution of P * A = transpose(S), avoiding an explicit inverse.
  pairing_2 = matrix(coefficient_ring, n, n, [
    integrate(basis[j] * basis[k]) for j in 1:n for k in 1:n
  ])
  pairing = fast_mode ? _specialize_matrix_at_origin(G, pairing_2) : pairing_2
  invariant_sums = zero_matrix(coefficient_ring, n, n)
  for degree in 0:(2 * valency(G))
    for beta in _effective_classes_with_functional_value(
        GKM_second_homology(G), first_chern_class(G), degree)
      show_progress && println("Computing quantum product in curve class $beta:")
      if iszero(beta)
        for j in 1:n, k in j:n
          invariant = integrate(class * basis[j] * basis[k])
          invariant_sums[j, k] += invariant
          j == k || (invariant_sums[k, j] += invariant)
        end
        continue
      end
      symmetric_indices = [(j, k) for j in 1:n for k in j:n]
      # marked_insertions = [
      #   ev(1, class) * ev(2, basis[j]) * ev(3, basis[k])
      #   for (j, k) in symmetric_indices
      # ]
      # invariants = gromov_witten(
      #   G, beta, 3, marked_insertions; show_bar=show_progress, fast_mode,
      # )
      marked_insertions = [
        GKMClass[class, basis[j], basis[k]]
        for (j, k) in symmetric_indices
      ]
      invariants = gromov_witten_nomarks(
        G, beta, marked_insertions; show_bar=show_progress, fast_mode,
      )
      for ((j, k), invariant) in zip(symmetric_indices, invariants)
        invariant_sums[j, k] += invariant
        j == k || (invariant_sums[k, j] += invariant)
      end
    end
  end
  result = _MP_inv(pairing) * invariant_sums
  fast_mode ? _specialize_matrix_at_origin(G, result) : result
end

function _specialize_matrix_at_origin(G::AbstractGKMGraph, M)
  origin = zeros(Int, rank_torus(G))
  matrix(QQ, nrows(M), ncols(M), [begin
    numerator_value = evaluate(numerator(M[i, j]), origin)
    denominator_value = evaluate(denominator(M[i, j]), origin)
    iszero(denominator_value) && throw(ArgumentError(
      "the multiplication matrix has no entrywise non-equivariant limit " *
      "in the chosen basis; use a basis that remains a basis at t=0",
    ))
    QQ(numerator_value) / QQ(denominator_value)
  end for i in 1:nrows(M) for j in 1:ncols(M)])
end
@doc raw"""
    quantum_product_at_q1(G::AbstractGKMGraph, class;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the matrix of equivariant quantum multiplication by `class` after
setting ``q=1``. The matrix represents

```math
  a \longmapsto (\mathrm{class} \ast a)|_{q=1}
```

in the chosen basis of ``H_T^*(X;\mathbb{Q})``, where `G` is the GKM graph of
``X``. The class can be supplied as a [`GKMClass`](@ref) or as a vector
containing one fixed-point restriction per vertex.

By default, the matrix is expressed in the basis
`[point_class(G, v) for v in 1:num_vertices(G)]`. Pass a collection of
`GKMClass` objects through `basis` to use a different cohomology basis. The
``j``-th column is the product of `class` with the ``j``-th basis element.

The graph must be strictly nef. Otherwise infinitely many effective curve
classes may contribute, so specialization to ``q=1`` need not be defined.
Set `show_progress=true` to display the underlying Gromov--Witten calculations.
Set `fast_mode=true` to use the faster graph-enumeration mode in positive
curve degrees and return the matrix after setting all equivariant parameters
to zero. The chosen basis must remain a basis under this specialization.

# Example

```jldoctest quantum_product_at_q1
julia> P1 = projective_space(GKMGraph, 1);

julia> quantum_product_at_q1(P1, point_class(P1, 1))
[(t1^2 - 2*t1*t2 + t2^2 + 1)//(t1 - t2)    1//(t1 - t2)]
[                         -1//(t1 - t2)   -1//(t1 - t2)]

julia> t1, t2 = gens_coeffRing(P1);

julia> quantum_product_at_q1(P1, [t1, t2])
[(t1^2 - t1*t2 + 1)//(t1 - t2)                    1//(t1 - t2)]
[                -1//(t1 - t2)   (t1*t2 - t2^2 - 1)//(t1 - t2)]
```
"""
function quantum_product_at_q1(G::AbstractGKMGraph, class;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  if class isa AbstractVector
    isempty(class) && throw(ArgumentError("class must not be empty"))
    class = localized_class(G, class)
  elseif class isa GKMClass
    class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  else
    throw(ArgumentError("class must be a GKMClass or a vector of fixed-point restrictions"))
  end
  _multiplication_matrix_at_q1(G, class; basis, show_progress, fast_mode)
end


@doc raw"""
    c1_at_q1(G::AbstractGKMGraph;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the matrix of equivariant quantum multiplication by ``c_1(T_X)`` at
``q=1``. By default, the point classes form the basis; pass a collection of
`GKMClass` objects through `basis` to use another cohomology basis.

The graph must be strictly nef so that only finitely many effective curve
classes contribute. Set `show_progress=true` to show the underlying
Gromov--Witten calculations. Set `fast_mode=true` to use the faster
graph-enumeration mode and return the matrix after setting all equivariant
parameters to zero. The chosen basis must have a well-defined non-equivariant
specialization.

# Example
```jldoctest c1_at_q1
julia> P1 = projective_space(GKMGraph, 1);

julia> c1_at_q1(P1)
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]
```
"""
function c1_at_q1(G::AbstractGKMGraph;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  quantum_product_at_q1(
    G, first_chern_class(G); basis, show_progress, fast_mode,
  )
end

@doc raw"""
    conjecture_O_eigenvalues(G::AbstractGKMGraph; printData::Bool=true)

Return the eigenvalues of quantum multiplication by ``c_1^T(TX)`` at ``q=1``
and at the non-equivariant limit ``t=0``.

The characteristic polynomial is computed before taking the limit, so this is
well-defined even when the standard-basis matrix contains rational functions
in the equivariant parameters. The graph must be compact and strictly nef.
When `printData=true`, the specialized characteristic polynomial is printed.

# Example
```jldoctest conjecture_O_eigenvalues
julia> P1 = projective_space(GKMGraph, 1);

julia> conjecture_O_eigenvalues(P1)
Characteristic poly of c1(TX)* at q=1, t=0:
x^2 - 4
2-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a1: -2.00000}
```
"""
function conjecture_O_eigenvalues(G::AbstractGKMGraph; printData::Bool=true)
  @req is_compact(G) "the graph must be compact to take the non-equivariant limit"
  chi = characteristic_polynomial(c1_at_q1(G))
  chi0 = polynomial(QQ, [0])
  z = zeros(Int, rank_torus(G))
  for i in 0:(length(chi) - 1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  printData && println("Characteristic poly of c1(TX)* at q=1, t=0:\n$chi0")
  roots(QQBar, chi0)
end

@doc raw"""
    twisted_c1_matrix(V::AbstractGKMVectorBundle, beta::CurveClass;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the ``q^\beta`` part of the matrix of the `V`-twisted equivariant
quantum product with ``c_1(T_X)-c_1(V)`` on ``X=\operatorname{baseof}(V)``.
The point classes are used by default; pass `basis` to choose another basis.

The twisting class is [`reduced_virtual_zero_section`](@ref). The resulting
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
  coefficient_ring = parent(first(localize(class).restrictions))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the base graph",
  ))
  pairing_2 = matrix(coefficient_ring, n, n, [
    integrate(euler_of_V * basis[j] * basis[k]) for j in 1:n for k in 1:n
  ])
  pairing = fast_mode ? _specialize_matrix_at_origin(G, pairing_2) : pairing_2
  invariant_sums = zero_matrix(coefficient_ring, n, n)
  if iszero(beta)
    for j in 1:n, k in j:n
      invariant = integrate(euler_of_V * class * basis[j] * basis[k])
      invariant_sums[j, k] += invariant
      j == k || (invariant_sums[k, j] += invariant)
    end
  else
    twist = virtual_zero_section(V)
    class_products = [GKMClass[class, basis[j]] for j in 1:n for k in 1:n]
    symmetric_indices = [(j, k) for j in 1:n for k in j:n]
    marked_insertions = [
      ev(1, class) * ev(2, basis[j]) * ev(3, basis[k]) * twist
      for (j, k) in symmetric_indices
    ]
    invariants = gromov_witten(
      G, beta, 3, marked_insertions; show_bar=show_progress, fast_mode,
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
  end
  # result = transpose(solve(transpose(pairing), invariant_sums))
  result = _MP_inv(pairing) * invariant_sums
  fast_mode ? _specialize_matrix_at_origin(G, result) : result
end

function _test_twisting_positivity_constraint(V::AbstractGKMVectorBundle)
  class = first_chern_class(baseof(V)) - first_chern_class(V)
  all(e -> begin d = integrate(class, e); denominator(d) == 1 && is_constant(numerator(d)) && constant_coefficient(numerator(d)) > 0 end, edges(baseof(V)))
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
      show_progress && println("Computing twisted c1 matrix in curve class $beta:")
      result += twisted_c1_matrix(V, beta; basis, show_progress, fast_mode)
    end
  end
  return result
  # chi = characteristic_polynomial(result)
  # chi0 = polynomial(QQ, [0])
  # z = zeros(Int, rank_torus(G))
  # for i in 0:(length(chi) - 1)
  #   set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  # end
  # (roots(QQBar, chi0), chi0, result)
end

function _MP_inv(A)
  m, n = nrows(A), ncols(A)
  coefficient_ring = base_ring(A)
  r, reduced = rref(A)
  iszero(r) && return zero_matrix(coefficient_ring, n, m)

  # Choose the pivot columns to obtain a rank factorization A = C * F,
  # with C of full column rank and F of full row rank.
  pivots = Vector{Int}(undef, r)
  for i in 1:r
    pivot = findfirst(j -> !iszero(reduced[i, j]), 1:n)
    isnothing(pivot) && error("failed to find a pivot in a nonzero RREF row")
    pivots[i] = pivot
  end
  C = A[:, pivots]
  F = transpose(solve(transpose(C), transpose(A)))

  # For A = C*F, the Moore--Penrose inverse is F^+*C^+.
  Ft = transpose(F)
  Ct = transpose(C)
  F_plus = solve(F * Ft, Ft)
  C_plus = transpose(solve(transpose(Ct * C), C))
  F_plus * C_plus
end
