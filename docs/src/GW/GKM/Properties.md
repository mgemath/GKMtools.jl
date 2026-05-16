# Properties of GKM graphs

These are some of the main properties of GKM graphs.

## General Properties
The following properties make sense for any abstract GKM graph.
```@docs
valency
rank_torus
is_compact
is2_indep
is3_indep
```

## Symplectic / almost complex properties
The following properties are well-defined as soon as there is a well-defined first Chern class.
This holds for GKM graphs coming from $T$-compatibly almost complex manifolds and from $T$-compatibly symplectic manifolds.
In the latter case, all choices of almost complex structures compatible with the symplectic form yield the same
first Chern class, making the properties below intrinsic to the symplectic structure.
```@docs
fano_index
pseudo_index
is_strictly_nef
```

## Compact Hamiltonian properties
The following properties have natural geometric meaning for GKM graphs of compact Hamiltonian GKM spaces.
That being said, they can be defined and checked for any abstract GKM graph, although one needs to be careful about their
geometric interpretation outside of the Hamiltonian setting.
```@docs
betti_numbers
index_periodic_betti
QH_ss_check_GLLXBR
```

## Index properties

```@docs
is_generic
xi_index
is_index_increasing
is_weakly_index_increasing
generic_xi_representatives
index_increasing_xi_representatives
weakly_index_increasing_xi_representatives
admits_index_increasing_xi
admits_weakly_index_increasing_xi
```