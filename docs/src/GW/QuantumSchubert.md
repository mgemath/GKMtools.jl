# Quantum Schubert calculus

For homogeneous spaces, this backend computes **ordinary small quantum**
Schubert structure constants from root data, without stable-map localization.
It uses equivariant quantum Chevalley and associativity internally, evaluated
at exact positive integer values of the simple roots. Only coefficients of
equivariant degree zero are returned; intermediate equivariant values must
not be interpreted as ordinary coefficients.

```julia
using Oscar, GKMtools

# E7/P2: no GKM graph or enumeration of the ambient Weyl group needed.
Q = quantum_schubert_context(root_system(:E, 7), [1, 3, 4, 5, 6, 7])
length(Q.representatives)  # 576
Q.novikov_degrees         # [14], in complex codimension

# D is the divisor associated with omitted simple root 2.
u = findfirst(==(3), Q.codimensions)
quantum_chevalley_product(Q, 2, u)
quantum_schubert_product(Q, u, u)
quantum_schubert_product(Q, u, u; max_degree=1)

# Explicitly request every product (potentially expensive).
# table = quantum_schubert_table(Q)
```

To preserve the vertex numbering of an existing graph, use
`Q = quantum_schubert_context(E7P2)`. The root-system constructor instead
sorts representatives by length and reduced word. Do not mix these index
orders. Inspect `Q.representatives` or pass a representative string.

The output uses positive geometric Schubert classes. With the current graph
axial-weight convention, the ordinary class represented by a Billey
`schubert_basis(G, v)` needs the factor `(-1)^Q.codimensions[v]` for this
normalization. This matters when comparing localization invariants (including
the third insertion / Poincaré dual), especially in odd dimension or odd
anticanonical degree. The graph constructor preserves labels, not this sign
of the `GKMClass` lift.

Each product is a sparse dictionary with keys `(w, d)`: vertex `w`
denotes the output Schubert class and tuple `d` its Novikov exponent.
Values are exact Oscar rational integers. For higher Picard rank, coordinates
are ordered by `Q.omitted_roots`, the increasing list of omitted simple roots.
These tuples are Schubert-curve coordinates, not arbitrary `CurveClass`
coordinates from a GKM graph.

For example, on projective two-space:

```julia
Q = quantum_schubert_context(root_system(:A, 2), [2])
quantum_schubert_product(Q, 3, 3) == Dict((2, (1,)) => QQ(1)) # h² ⋆ h² = qh
quantum_schubert_coefficient(Q, 3, 3, 2, 1) == 1
```

All degrees allowed by grading are included by default. The optional
`max_degree` limits returned Novikov exponents componentwise; it does not
discard intermediate equivariant coefficients required for reconstruction.
For E7/P2, dimension 42 and index 14 give a global bound of degree 6 for
ordinary basis products. The context caches coefficients lazily. A complete
table can still require substantial time and memory; no E7/P2 full-table
performance guarantee is made.

The module is independent of the existing localization-based
`quantum_product` API. It does not change multiplication of `GKMClass`
objects, which remains classical, or automatically replace generic GKM
quantum calculations. It applies to complete homogeneous spaces G/P, not to
arbitrary Schubert subvarieties, twisted theories, or higher-genus invariants.

## Algorithm

Sparse divisor edges use positive roots outside the Levi, their coroot
coordinates, and parabolic length conditions. The classical target must
itself be minimal; the quantum target is reduced to its minimal coset
representative. Coroots, rather than root coefficients, are necessary in
non-simply-laced types.

For reconstruction the ample divisor is the sum of the Schubert divisors.
Its equivariant diagonal is obtained from Billey's reduced-word formula.
A sufficiently large integer radix gives distinct diagonal evaluations at
all vertices and positive values on positive roots. Associativity then gives
a recurrence with nonzero rational denominators. Classical steps increase
equivariant degree; quantum steps decrease the effective multidegree.
Diagonal coefficients are solved from the unit equation as affine expressions
in one unknown. Commutativity and memoization avoid repeated calculations.

The mathematical source is Mihalcea,
[Equivariant quantum cohomology of homogeneous spaces](https://arxiv.org/abs/math/0501213),
especially Corollary 6.5 and Sections 7–8. The exact specialization is an
implementation choice: degree-zero equivariant polynomials are constants, so
their evaluated values give the ordinary invariants without interpolation.

## API

```@docs
QuantumSchubertContext
quantum_schubert_context
quantum_chevalley_product
quantum_schubert_coefficient
quantum_schubert_product
quantum_schubert_table
serialize_quantum_schubert_products
quantum_schubert_matrix
```
