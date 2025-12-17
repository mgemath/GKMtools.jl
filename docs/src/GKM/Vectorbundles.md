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