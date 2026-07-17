# Standard Constructions

In this section we present some functions that allow the construction of famous GKM graphs.

## Generalized Flags Varieties

```@docs
generalized_gkm_flag
flag_variety
grassmannian
projective_space
```

## Generalized Schubert Varieties
```@docs
get_bruhat_order_of_generalized_flag
generalized_gkm_schubert
bott_samelson
```
## Schubert Classes

By a _Schubert class_ we simply mean the output of [`poincare_dual`](@ref) applied to a Schubert variety, considered as
subvariety of another Schubert variety (which could be a homogeneous space).

!!! note
    When a Schubert variety is not smooth, its GKM graph is not necessarily regular.
    It follows from the definition of [`poincare_dual`](@ref) that the returned Schubert class
    may have different degrees at different vertices if the Schubert subvariety or the 
    containing Schubert variety is not smooth.
    That is, the localization of the Schubert class at each vertex is a homogeneous polynomial,
    but the degrees of these polynomials could vary as the vertex varies.


```@docs
schubert_class
schubert_classes
```

## Toric varieties

```@docs
gkm_graph_of_toric
```