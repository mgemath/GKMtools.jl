##############################################
###   Potential optimizations in future:  ###
##############################################
# 1. Use the nomarks feature
# 2. Use fast mode. There are two possible ways:
#     - A. Use the same random t1,...,tr throughout.
#          This could cause random eigenvalues. Still useful for atoms?
#     - B. Work with bases for equivariant cohomology that are actual bases,
#          i.e., without localizing the coefficient ring.
#          Then we can restrict to using GW invariants that are in QQ, so the
#          current fastmode works.

#     - Note: I think fast mode and nomarks are not compatible yet.
#             Nomarks works whenever we have no psi and lambda classes,
#             which should be easy enough to check on an Equivariant Class.
##############################################


@doc raw"""
    twisted_c1_matrix(V::GKM_vector_bundle, beta; show_progress::Bool=false)

Return the $q^\beta$ part of the matrix given by the equivariant quantum product on $X$ with $c_1(T_X) - c_1(V)$, twisted by the convex vector bundle $V$.
Here, $X$ is the base of the vector bundle $V$.

!!! note
    See [MR1685628; Section 2.1](@cite) for the definition of the twisted quantum product, and see also [`reduced_virtual_zero_section`](@ref) for an implementation of general twisted quantum products.
    The twisted quantum cohomology on $X$ admits a ring homomorphism to the ordinary quantum cohomology of the smooth zero-locus $Y\subset X$ of a section of $V$.
    Since we really care about the quantum product by $c_1(T_Y)$ on $Y$, we need some class on $X$ that restricts to $c_1(T_Y))$.
    This is precisely $c_1(T_X) - c_1(V)$.

# Input
- `V::GKM_vector_bundle`: A _convex_ vector bundle.
- `beta`: A curve class on the GKM graph `baseof(V)`, the base of `V`.
- `show_progress::Bool` (optional): If set to true, the progress bars of [`gromov_witten`](@ref) will be shown.

# Output
The matrix of the $q^\beta$ part of the quantum product by $c_1(T_X) - c_1(V)$, twisted by $V$.

!!! note
    The returned matrix is expressed in the basis $(1, 0, \ldots, 0), (0, 1, 0,\ldots, 0), \ldots, (0,\ldots,0,1)$ of $H_T^*(X;\mathbb{Q})$ localized at the 
    fraction field of the coefficient ring. These classes do not represent classes in $H_T^*(X;\mathbb{Q})$ without localizing the coefficient ring,
    so in particular the output will consist of rational functions even when `V` comes from a vector bundle on a projective GKM space.

# Example
As an example, let us see the case $V=\mathcal{O}_{\mathbb{P}^2}(1)$ in degrees up to $2$.
```@jldoctest
julia> V = vector_bundle_O(2, [1]);

julia> beta = curve_class(baseof(V), "1", "2");

julia> twisted_c1_matrix(V, 0*beta)
[2*t1 - t2 - t3 - t4              0               0]
[                  0   t2 - t3 - t4               0]
[                  0              0   -t2 + t3 - t4]

julia> twisted_c1_matrix(V, 1*beta)
[               (2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)                  (2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)                  (2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)]
[ (2*t1 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)    (2*t1 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)    (2*t1 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)]
[(-2*t1 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)   (-2*t1 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)   (-2*t1 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)]

julia> twisted_c1_matrix(V, 2*beta)
[0   0   0]
[0   0   0]
[0   0   0]
```
Note that the returned matrices contain all equivariant parameters for `V`, including those that act only in the fibre direction.
"""
function twisted_c1_matrix(V::GKM_vector_bundle, beta; show_progress::Bool=false)

  G = baseof(V)
  nv = n_vertices(G.g)

  M = zero_matrix(G.equivariantCohomology.coeffRingLocalized, nv, nv)
  g = gens(G.equivariantCohomology.cohomRingLocalized)

  P = reduced_virtual_zero_section(V)
  class = chern_class(G, 1) - chern_class(V, 1)

  for i in 1:nv
    cg = quantum_product(G, beta, class, g[i]; show_progress=show_progress, useStructureConstants=false, twist_class=P)
    M[i,:] += [cg[j] for j in 1:nv]
  end

  return M
end


@doc raw"""
    twisted_c1_matrix_at_q1(V::GKM_vector_bundle; show_progress::Bool=false)

Return the matrix of the equivariant twisted quantum product by $c_1(T_X) - c_1(V)$ on $X$, twisted by the _convex_ vector bundle $V$ where $X$ is the base of $V$.
This sums the $q^\beta$ terms produced by [`twisted_c1_matrix`](@ref), setting $q=1$ and summing over all relevant curve classes $\beta$.

The output also contains the characteristic polynomial of the resulting matrix and its eigenvalues after setting all equivariant parameters to zero.

!!! note
    This function requires that the second cohomology class $c_1(T_X) - c_1(V)$ is positive on all $T$-stable copies of $\mathbb{P}^1$ inside $X$.
    Otherwise, there is no guarantee that only finitely many terms contribute.

# Input
- `V::GKM_vector_bundle`: A _convex_ vector bundle.

# Output
The output is a tuple `(roots, chi0, M)`, where:
- `M` is the described matrix.
- `chi0` is the characteristic polynomial with all equivariant parameters set to zero.
- `roots` are the complex roots of `chi0`, i.e. the eigenvalues in the non-equivariant limit of `M`.

!!! note
    See [`twisted_c1_matrix`](@ref) for the choice of basis in which `M` is returned.
    This basis is not a basis over $H_T^*$, but only over $\text{Frac}(H_T^*)$, so
    the entries of `M` are rational functions in the equivariant parameters.
    However, since the characteristic polynomial $\chi_M$ of $M$ is independent of the choice of basis for `M`, its coefficients are polynomials in the equivariant parameters.
    In particular, the non-equivariant limit `chi0` of $\chi_M$ is well-defined when `V` is realized as convex vector bundle over a projective GKM space.

# Example 1

Let us first continue the example of $V=\mathcal{O}_{\mathbb{P}^2}(1)$ from [`twisted_c1_matrix`](@ref).

```@jldoctest
julia> V = vector_bundle_O(2, [1]);

julia> (roots, chi0, M) = twisted_c1_matrix_at_q1(V);

julia> roots
3-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a1: 0}
 {a1: -2.00000}

julia> chi0
x^3 - 4*x

julia> M
[(2*t1^3 - 3*t1^2*t2 - 3*t1^2*t3 - t1^2*t4 + t1*t2^2 + 4*t1*t2*t3 + t1*t2*t4 + t1*t3^2 + t1*t3*t4 - t2^2*t3 - t2*t3^2 - t2*t3*t4 + 2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)                                                                                                                          (2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)                                                                                                                           (2*t4)//(t1^2 - t1*t2 - t1*t3 + t2*t3)]
[                                                                                                                   (2*t1 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)   (t1*t2^2 - 2*t1*t2*t3 - t1*t2*t4 + t1*t3^2 + t1*t3*t4 + 2*t1 - t2^3 + 2*t2^2*t3 + t2^2*t4 - t2*t3^2 - t2*t3*t4 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)                                                                                                             (2*t1 - 2*t2 - 2*t4)//(t1*t2 - t1*t3 - t2^2 + t2*t3)]
[                                                                                                                  (-2*t1 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)                                                                                                           (-2*t1 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)   (-t1*t2^2 + 2*t1*t2*t3 - t1*t2*t4 - t1*t3^2 + t1*t3*t4 - 2*t1 + t2^2*t3 - 2*t2*t3^2 + t2*t3*t4 + t3^3 - t3^2*t4 + 2*t3 + 2*t4)//(t1*t2 - t1*t3 - t2*t3 + t3^2)]
```

# Example 2
Next, let us see $\mathcal{O}(1)$, $\mathcal{O}(1)\oplus\mathcal{O}(1)$, and $\mathcal{O}(2)$ on $\mathbb{P}^3$.
We begin with $\mathcal{O}_{\mathbb{P}^3}(1)$.

```jldoctest
julia> V = vector_bundle_O(3, [1]);

julia> (roots, chi0, M) = twisted_c1_matrix_at_q1(V);

julia> chi0
x^4 - 27*x

julia> roots
4-element Vector{QQBarFieldElem}:
 {a1: 3.00000}
 {a1: 0}
 {a2: -1.50000 + 2.59808*im}
 {a2: -1.50000 - 2.59808*im}
```

Next, let us see $\mathcal{O}_{\mathbb{P}^3}(1)\oplus\mathcal{O}_{\mathbb{P}^3}(1)$.

```jldoctest
julia> V = vector_bundle_O(3, [1, 1]);

julia> (roots, chi0, M) = twisted_c1_matrix_at_q1(V);

julia> chi0
x^4 - 4*x^2

julia> roots
4-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a1: 0}
 {a1: 0}
 {a1: -2.00000}
```

Finally, let us see $\mathcal{O}_{\mathbb{P}^3}(2)$.

```jldoctest
julia> V = vector_bundle_O(3, [2]);

julia> (roots, chi0, M) = twisted_c1_matrix_at_q1(V);

julia> chi0
x^4 - 16*x^2

julia> roots
4-element Vector{QQBarFieldElem}:
 {a1: 4.00000}
 {a1: 0}
 {a1: 0}
 {a1: -4.00000}
```
"""
function twisted_c1_matrix_at_q1(V::GKM_vector_bundle; show_progress::Bool=false)

  G = baseof(V)
  nv = n_vertices(G.g)
  @req _test_twisting_positivity_constraint(V) "c1(G)-c1(V) needs to be strictly positive on all edges"
  
  max_functional_number = 2*valency(G)
  M = zero_matrix(G.equivariantCohomology.coeffRingLocalized, nv, nv)

  for c in 0:max_functional_number
    for beta in _effective_classes_with_functional_value(GKM_second_homology(G), chern_class(G, 1) - chern_class(V, 1), c)

      show_progress && println("Computing twisted c1 matrix in curve class $beta:")
      M += twisted_c1_matrix(V, beta; show_progress=show_progress)
    end
  end

  chi = characteristic_polynomial(M)
  chi0 = polynomial(QQ, [0])
  z = repeat([0], rank_torus(G))
  for i in 0:(length(chi)-1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  return (roots(QQBar, chi0), chi0, M)
end

function _test_twisting_positivity_constraint(V::GKM_vector_bundle)::Bool
  G = baseof(V)
  dc = chern_class(G, 1) - chern_class(V, 1)
  for e in edges(G.g)
    de_frac_type = integrate(dc, G, e)
    @req denominator(de_frac_type) == 1 "difference of first chern classes of G and v is not a GKM glass"
    @req is_constant(numerator(de_frac_type)) "degree over e is not an integer! This should not happen."
    de = constant_coefficient(numerator(de_frac_type))
    de <= 0 && return false
  end
  return true
end
