# Examples

The following pages provide the main tools to construct GKM graphs, with some additional tools (e.g., vector bundles). In particular, we provided:

## Generalized flag varieties
In particular, we provide two different ways to construct those varieties:
* [`generalized_gkm_flag`](@ref): the most complete one, it takes as input a root system,
* [`flag_variety`](@ref): it returns a flag varieties, vertices are indexed by flags.

If the variety is constructed using [`generalized_gkm_flag`](@ref), we provide additional tools:
* The constructor of tautological vector bundles, [`tautological_bd`](@ref).
* The constructor of the cohomology ring as a Oscar graded ring, with Schubert basis, [`schubert_cohomology_ring`](@ref).

## Toric varieties
The main source of orbifold GKM graphs. We have two type of toric varieties.

### Smooth toric varieties
It supports Oscar toric varieties and line bundles. Those can be converted to GKM smooth graphs and vector bundles. Use [`gkm_graph_of_toric`](@ref).

### Orbifold toric varieties
The main examples are weighted projective spaces. Moreover, the tautological line bundle $\mathcal{O}(1)$ can be constructed.
In order to obtain an orbifold toric variety from a Oscar toric variety, use [`gkm_graph_of_orbifold_toric`](@ref).

Oscar supports weighted projective space, but withour the orbifold structure. That is, using Oscar:
```julia
julia> using Oscar

julia> X = weighted_projective_space(NormalToricVariety, [1, 2, 4]);

julia> Y = weighted_projective_space(NormalToricVariety, [1, 1, 2]);
```
The objects $X$ and $Y$ are the same, thus calling
```julia
julia> gkm_graph_of_orbifold_toric(X);

julia> gkm_graph_of_orbifold_toric(Y);
```
produces the same object, that is $\mathbb{P}(1, 1, 2)$. In order to obtain $\mathbb{P}(1, 2, 4)$, use:
```julia
julia> fan = stacky_weighted_projective_space_fan([1, 2, 4]);

julia> P124 = gkm_graph_of_orbifold_toric(fan);
```

## Other example
We provide the constructors [`gkm_2d`](@ref), [`gkm_3d_twisted_flag`](@ref) and [`gkm_3d_positive_non_toric`](@ref) for interesting examples in low dimension.