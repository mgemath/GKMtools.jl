```@meta
DocTestSetup = :(using Oscar, GKMtools)
```

# Generalized Flag Varieties

```@docs
generalized_gkm_flag
flag_variety
grassmannian
projective_space
```

## Generalized Schubert Varieties

[`generalized_gkm_schubert`](@ref) constructs the induced GKM subgraph of a
**smooth** Schubert variety. With `descending=false`, it uses the opposite
Schubert variety instead. The constructor also computes smooth GKM homology
data, so it is not a general constructor for singular Schubert moment graphs.

There are two distinct limitations for singular varieties:

- If the moment graph has nonconstant vertex valency, the Schubert variety is
  **not rationally smooth**. The constructor throws an `ArgumentError` explaining
  this obstruction before attempting the smooth GKM construction.
- A **rationally smooth but singular** Schubert variety has constant vertex
  valency, which alone is insufficient for smoothness. Construction can still
  fail, for example when Chern weights are not proportional to an edge weight.
  The `ArgumentError` explains this distinction and retains the underlying error.

A successful return does not certify smoothness. Establish smoothness separately
before interpreting the graph's invariants geometrically. To compute Schubert
classes in the smooth ambient flag variety, including classes indexed by singular
Schubert varieties, use [`schubert_basis`](@ref).

For example, the following calls illustrate the two failure cases:

```jldoctest
julia> generalized_gkm_schubert(root_system(:A, 3), "s2*s1*s3*s2");
ERROR: ArgumentError: The Schubert variety indexed by s2*s1*s3*s2 is not rationally smooth: its moment graph does not have constant valency. generalized_gkm_schubert supports smooth Schubert varieties only.

julia> generalized_gkm_schubert(root_system(:B, 2), "s2*s1*s2");
ERROR: ArgumentError: Cannot construct the smooth GKM graph of the Schubert variety indexed by s2*s1*s2. Its moment graph has constant valency, but rational smoothness does not imply smoothness. Singular Schubert varieties are not supported by generalized_gkm_schubert. Underlying error: ArgumentError: Chern weights are not proportional to the edge weight
```

```@docs
generalized_gkm_schubert
bott_samelson
```

## Schubert Classes
```@docs
schubert_basis
billey_schubert_basis
schubert_class
schubert_classes
```

## Cohomology Ring
```@docs
schubert_cohomology_ring
```

## Vector bundles on G/P

Let $R$ be a root system, and $S$ a (possibly empty) set of simple roots of $R$ defining a parabolic subgroup $P\subseteq G$, where $G$ is the simply connected Lie group defined by $R$. See [`generalized_gkm_flag`](@ref) for the construction of the GKM graph of $G/P$.

Following [MR89473](@cite), there is a natural correspondence between $G$-equivariant bundles over $G/P$ and representations $(\rho, V)$ of $P$. Given a representation $\rho\colon P \rightarrow GL(V)$, the relative vector bundle is given by
```math
  G\times_P V, \text{ where  }(g, v) \sim (gp, \rho(p)^{-1}(v)).
```
Irreducible finite dimensional representations of $P$ are in one-to-one correspondence with $L$-dominant weights $\lambda$, that is with weights of $G$ such that $\langle \lambda, \alpha^{\vee}\rangle \ge 0$ for all $\alpha \in S$.

We provide a function that constructs the vector bundle corresponding to any $\lambda$. We define weights using Oscar's [functions](https://docs.oscar-system.org/stable/LieTheory/weight_lattices/).

```@docs
rank_of_bd
tautological_bd
```
