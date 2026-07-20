# Vector Bundles

Theory of GKM vector bundles, see [GSZ12](@cite).

## Construction

```@docs
line_bundle
vector_bundle
orbifold_vector_bundle
orbifold_line_bundle
tangent_bundle
fiber_weight
fiber_representation
```

## GKM graphs from vector bundles

```@docs
baseof
```

## Linear algebra

```@docs
direct_sum
tensor_product
wedge_product
sym_product
```

## Connections

```@docs
build_vector_bundle_connection
```

## Important Examples

```@docs
line_bundle_O
weighted_projective_line_bundle
tautological_and_univ_bd
```

### Generalized GKM flags

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