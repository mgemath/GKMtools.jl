# Vector Bundles

Theory of GKM vector bundles, see [GSZ12](@cite).

## Construction

```@docs
line_bundle
vector_bundle
tangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
cotangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
projectivization
```

## Linear algebra

```@docs
rank(::GKMtools.GKM_vector_bundle)
direct_sum
dual(::GKMtools.GKM_vector_bundle)
*(::GKMtools.GKM_vector_bundle, ::GKMtools.GKM_vector_bundle)
^(::GKMtools.GKM_vector_bundle, ::Number)
wedge_product
sym_product
```

## Connections

```@docs
get_connection(::GKMtools.GKM_vector_bundle)
get_any_connection(::GKMtools.GKM_vector_bundle)
```

## Important Examples

```@docs
vector_bundle_O(::Int64, ::Vector{Int64})
univ_quotient_bd
tautological_bd
```