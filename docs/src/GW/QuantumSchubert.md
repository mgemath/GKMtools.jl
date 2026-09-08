# Quantum Schubert calculus

For homogeneous spaces, this backend computes **ordinary small quantum**
Schubert structure constants from root data, without stable-map localization.
It uses equivariant quantum Chevalley and associativity internally, evaluated
at admissible exact rational or finite-field values of the simple roots. Only coefficients of
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
Values are exact Oscar rational integers in an exact context, or finite-field
elements in a modular context. For higher Picard rank, coordinates
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
Diagonal coefficients are solved from a scalar unit recurrence using
precomputed classical diagonal restrictions. Commutativity and memoization avoid repeated calculations.

The mathematical source is Mihalcea,
[Equivariant quantum cohomology of homogeneous spaces](https://arxiv.org/abs/math/0501213),
especially Corollary 6.5 and Sections 7–8. The exact specialization is an
implementation choice: degree-zero equivariant polynomials are constants, so
their evaluated values give the ordinary invariants without interpolation.


## Large spaces and computation modulo a prime

Use `p` to perform coefficient reconstruction directly in a finite field,
without first computing large rational numbers:

```julia
p = 1_000_000_007
Q = quantum_schubert_context(E7P5; p, max_cache_entries=200_000)
serialize_quantum_schubert_products("$namefile-$v-mod$p.jls", Q, v; show_bar=true)

using Serialization
products = deserialize("$namefile-$v-mod$p.jls")
M = quantum_schubert_matrix(products; p, at_q1=true)
```

The graph overload preserves your existing vertex indices. If a graph is not
otherwise needed, avoid constructing it:

```julia
Q = quantum_schubert_context(
    root_system(:E,7), [1,2,3,4,6,7];
    p=1_000_000_007, max_cache_entries=200_000,
)
# Choose v from Q.representatives / Q.codimensions in this ordering.
```

Omitting `p` retains exact rational arithmetic. With `p`, public coefficients
are finite-field elements. Serialization still stores only the product vector,
using `BigInt` representatives in `0:p-1`, with no metadata. In particular,
**the prime is not stored**: supply it again to `quantum_schubert_matrix`.
The context overload `quantum_schubert_matrix(Q,v)` uses `Q.prime`
automatically. Residues do not determine the original integer invariants;
this API does not perform Chinese-remainder or integer reconstruction.

The prime must exceed the number of Schubert classes. The constructor checks
primality and tries deterministic equivariant specializations, accepting one
only when all divisor diagonals are distinct and every positive-root value is
nonzero. It rejects the prime if no admissible specialization is found. A large
prime is preferable; these checks prevent invalid divisions modulo p.

The default `max_cache_entries=200_000` bounds each reconstruction memo,
including the coefficient, diagonal, and classical-restriction caches.
When full, a memo is cleared and subsequently repopulated; this trades some
recomputation for bounded retained entries. Setting it to zero disables
memoization. This is an **entry limit, not a total RAM limit**: root data,
Bruhat bitsets, the returned product vector, temporary allocations, and the
sizes of exact rational numbers still consume memory. A smaller limit may
increase computation time. The classical Bruhat index is built from cover
edges using bitsets and is omitted when its payload would exceed 64 MiB.

Classical coefficients with an input equal to the output are evaluated
directly by a backward Billey recursion with bounded memoization, rather than
expanded along upward Chevalley paths.

Diagonal reconstruction now uses a scalar unit recurrence and precomputed
classical diagonal restrictions, avoiding temporary polynomials. It routes
branches that cannot reach the diagonal pivot through the shared cache.
The associativity recurrence expands the higher-codimension input to reduce
unnecessary upward steps. Products only visit output classes in the required
codimension. Native prime-field elements avoid generic extension-field overhead. These changes
apply to both exact and modular contexts.

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
