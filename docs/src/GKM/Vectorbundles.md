# Vector Bundles

Theory of GKM vector bundles, see [GSZ12](@cite).

## Construction

```@docs
line_bundle
vector_bundle
tangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
cotangent_bd(::GKMtools.AbstractGKM_graph; ::Int64)
gkm_line_bundle_of_toric
gkm_vector_bundle_of_toric
```

## GKM graphs from vector bundles

```@docs
total_space(::GKMtools.GKM_vector_bundle)
projectivization
baseof
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
isvalid(::Dict{Tuple{Oscar.Edge, Int64}, Int64}, ::GKMtools.GKM_vector_bundle; ::Bool)
```

## Important Examples

```@docs
vector_bundle_O(::Int64, ::Vector{Int64})
tautological_and_univ_bd
```