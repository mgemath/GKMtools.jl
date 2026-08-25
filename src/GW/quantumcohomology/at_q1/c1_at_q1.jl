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
    show_progress::Bool=false)
  @req all(e -> chern_number(e, G) > 0, edges(G)) "the graph must be strictly nef"
  n = num_vertices(G)
  result = zero_matrix(parent(first(localize(class).restrictions)), n, n)
  basis = gens_cohomRing(G)
  for degree in 0:(2 * valency(G))
    for beta in _effective_classes_with_functional_value(
        GKM_second_homology(G), first_chern_class(G), degree)
      show_progress && println("Computing quantum product in curve class $beta:")
      for i in 1:n
        result[i, :] += (iszero(beta) ? localize(class * basis[i]).restrictions : gromov_witten(G, beta, 3, [ev(1, class) * ev(2, basis[i]) * ev(3, point_class(G, v)) for v in 1:n]; show_bar=show_progress))
      end
    end
  end
  result
end
@doc raw"""
    quantum_product_at_q1(G::AbstractGKMGraph, class; show_progress::Bool=false)

Return the matrix of equivariant quantum multiplication by `class` after
setting ``q=1``. The matrix represents

```math
  a \longmapsto (\mathrm{class} \ast a)|_{q=1}
```

in the [standard basis](#The-standard-basis) of ``H_T^*(X;\mathbb{Q})``, where
`G` is the GKM graph of ``X``. The class can be supplied as a [`GKMClass`](@ref)
or as a vector containing one fixed-point restriction per vertex.

The graph must be strictly nef. Otherwise infinitely many effective curve
classes may contribute, so specialization to ``q=1`` need not be defined.
Set `show_progress=true` to display the underlying Gromov--Witten calculations.

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
function quantum_product_at_q1(G::AbstractGKMGraph, class; show_progress::Bool=false)
  if class isa AbstractVector
    isempty(class) && throw(ArgumentError("class must not be empty"))
    class = localized_class(G, class)
  elseif class isa GKMClass
    class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  else
    throw(ArgumentError("class must be a GKMClass or a vector of fixed-point restrictions"))
  end
  _multiplication_matrix_at_q1(G, class; show_progress)
end


@doc raw"""
    c1_at_q1(G::AbstractGKMGraph; show_progress::Bool=false)

Return the matrix of equivariant quantum multiplication by ``c_1(T_X)`` at
``q=1`` in the [standard basis](#The-standard-basis).

The graph must be strictly nef so that only finitely many effective curve
classes contribute. Set `show_progress=true` to show the underlying
Gromov--Witten calculations.

# Example
```jldoctest c1_at_q1
julia> P1 = projective_space(GKMGraph, 1);

julia> c1_at_q1(P1)
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]
```
"""
function c1_at_q1(G::AbstractGKMGraph; show_progress::Bool=false)
  quantum_product_at_q1(G, first_chern_class(G); show_progress)
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
        show_progress::Bool=false)

Return the ``q^\beta`` part of the matrix of the `V`-twisted equivariant
quantum product with ``c_1(T_X)-c_1(V)`` on ``X=\operatorname{baseof}(V)``, in
the [standard basis](#The-standard-basis).

The twisting class is [`reduced_virtual_zero_section`](@ref). The resulting
quantum cohomology maps to that of the smooth zero locus ``Y\subset X`` of a
section of `V`, and ``c_1(T_X)-c_1(V)`` restricts to ``c_1(T_Y)``.

# Arguments
- `V`: a convex GKM vector bundle.
- `beta`: a curve class on `baseof(V)`; the zero class is allowed.
- `show_progress`: display progress for the Gromov--Witten computations.

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
    show_progress::Bool=false)
  G = baseof(V)
  _check_curve_class(G, beta)
  n = num_vertices(G)
  class = first_chern_class(G) - first_chern_class(V)
  result = zero_matrix(parent(first(localize(class).restrictions)), n, n)
  basis = gens_cohomRing(G)
  if iszero(beta)
    for i in 1:n
      result[i, :] += localize(class * basis[i]).restrictions
    end
    return result
  end
  twist = reduced_virtual_zero_section(V)
  for i in 1:n
    insertions = [twist * ev(1, class) * ev(2, basis[i]) * ev(3, point_class(G, v))
      for v in 1:n]
    result[i, :] += gromov_witten(G, beta, 3, insertions; show_bar=show_progress)
  end
  result
end

function _test_twisting_positivity_constraint(V::AbstractGKMVectorBundle)
  class = first_chern_class(baseof(V)) - first_chern_class(V)
  all(e -> begin d = integrate(class, e); denominator(d) == 1 && is_constant(numerator(d)) && constant_coefficient(numerator(d)) > 0 end, edges(baseof(V)))
end

@doc raw"""
    twisted_c1_matrix_at_q1(V::AbstractGKMVectorBundle;
        show_progress::Bool=false)

Return the matrix of the `V`-twisted equivariant quantum product with
``c_1(T_X)-c_1(V)`` at ``q=1``, together with its non-equivariant
characteristic polynomial and eigenvalues.

The function sums [`twisted_c1_matrix`](@ref) over all contributing effective
curve classes. It requires ``c_1(T_X)-c_1(V)`` to be strictly positive on every
invariant edge, ensuring that the sum at ``q=1`` is finite.

# Output
The result is `(roots, chi0, M)`, where `M` is the multiplication matrix,
`chi0` is its characteristic polynomial after setting the equivariant
parameters to zero, and `roots` are the roots of `chi0` in `QQBar`.

Although `M` can contain rational functions in the standard basis, its
characteristic polynomial is basis-independent, so its non-equivariant limit
is well-defined for a convex bundle over a projective GKM space.

# Example
```jldoctest twisted_c1_matrix_at_q1
julia> V = vector_bundle_O(1, [1]);

julia> roots, chi0, M = twisted_c1_matrix_at_q1(V);

julia> chi0
x^2 + x
```
"""
function twisted_c1_matrix_at_q1(V::AbstractGKMVectorBundle; show_progress::Bool=false)
  G = baseof(V)
  @req _test_twisting_positivity_constraint(V) "c1(G)-c1(V) must be strictly positive on every edge"
  class = first_chern_class(G) - first_chern_class(V)
  n = num_vertices(G)
  result = zero_matrix(parent(first(localize(class).restrictions)), n, n)
  for degree in 0:(2 * valency(G))
    for beta in _effective_classes_with_functional_value(GKM_second_homology(G), class, degree)
      show_progress && println("Computing twisted c1 matrix in curve class $beta:")
      result += twisted_c1_matrix(V, beta; show_progress)
    end
  end
  chi = characteristic_polynomial(result)
  chi0 = polynomial(QQ, [0])
  z = zeros(Int, rank_torus(G))
  for i in 0:(length(chi) - 1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  (roots(QQBar, chi0), chi0, result)
end
