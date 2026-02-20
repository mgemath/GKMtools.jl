# Vector Bundles

Theory of GKM vector bundles, see [GSZ12](@cite).

## Construction

```@docs
line_bundle
vector_bundle
tangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
cotangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
projectivization
baseof
gkm_line_bundle_of_toric
gkm_vector_bundle_of_toric
```

## Linear algebra

```@docs
rank(::GKMtools.GKM_vector_bundle)
direct_sum
+(::GKMtools.GKM_vector_bundle, ::GKMtools.GKM_vector_bundle)
dual(::GKMtools.GKM_vector_bundle)
*(::GKMtools.GKM_vector_bundle, ::GKMtools.GKM_vector_bundle)
^(::GKMtools.GKM_vector_bundle, ::Number)
wedge_product
sym_product
det(::GKMtools.GKM_vector_bundle)
```

## Connections

```@docs
get_connection(::GKMtools.GKM_vector_bundle)
get_any_connection(::GKMtools.GKM_vector_bundle)
```

## Important Examples

```@docs
vector_bundle_O(::Int64, ::Vector{Int64})
tautological_and_univ_bd
```

### Generalized GKM flags

Let $R$ be a root system, and $S$ a (possibly empty) set of simple roots of $R$ defining a parabolic subgroup $P\subseteq G$, where $G$ is the simply connected Lie group defined by $R$. See [`generalized_gkm_flag`](@ref) for the construction of the GKM graph of $G/P$.

Following [MR89473](@cite), there is a natural correspondence between $G$-equivariant bundles over $G/P$ and representations $(\rho, V)$ of $P$. Given a representation $\rho\colon P \rightarrow GL(V)$, the relative vector bundle is given by
```math
  G\times_P V, \text{ where  }(g, v) \sim (gp, \rho(p)^{-1}(v)).
```
Irreducible finite dimensional representations of $P$ are in one-to-one correspondence with $L$-dominant weights $\lambda$, that is with weights of $G$ such that $\langle \lambda, \alpha^{\vee}\rangle$ for all $\alpha \in S$.

We provide a function that constructs the vector bundle corresponding to any $\lambda$. We define weights using Oscar's [functions](https://docs.oscar-system.org/stable/LieTheory/weight_lattices/).

```@docs
rank_of_bd
tautological_bd
```