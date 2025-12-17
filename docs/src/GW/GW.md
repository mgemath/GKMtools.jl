
# Equivariant Gromov--Witten invariants of GKM graphs

One of the main features of this package is to calculate equivariant Gromov-Witten invariants of GKM spaces $X$ (see [Definition](../GKM/GKM.md#Definition)).

**Mathematical ingredients:**
The following mathematical ingredients are used in the package:
  - The localization formula from [LS17, Hirschi_2023](@cite) that expresses equivariant Gromov-Witten invariants in the form
    ```math
      GW^T_{X,g,m,\beta}(...) = \sum_{\overrightarrow{\Gamma}} GW^T_{\overrightarrow{\Gamma}}
    ```
    where the sum is over certain decorated graphs (trees in genus zero) $\overrightarrow{\Gamma}$ mapping into the GKM graph of $X$, such that the edge multiplicities sum to the curve class $\beta\in H_2(X;\mathbb{Z})$.
    Here, $GW^T_{\overrightarrow{\Gamma}}$ is a rational function in $\dim{T}$ many variables, i.e. an element of $\text{Frac}H_T^*(\text{pt};\mathbb{Q})$.

  - An efficient way of enumerating graphs (or trees in genus zero) with marked points and counting their automorphisms (inspired by [WROM86](@cite)).

  - A way of enumerating all combination of edges with multiplicities in the GKM graph that sum to a given curve class $\beta\in H_2(X;\mathbb{Z})$ (*see supporting paper*).

  - For positive genus, a table of Hodge integrals $\int_{\overline{\mathcal{M}}_{g,n}} \psi_1^{a_1}\cdots \psi_n^{a_n}\lambda_1^{l_1}\cdots\lambda_g^{l_g}$ up to the required genus, which we generate using [Yang2010, HodgeIntegralsSource](@cite).

Interesting aspects of the positive genus calculation will be discussed in a follow-up paper.

## Integrating on the moduli space

```@docs
gromov_witten
```

## Cohomology classes on the moduli space

```@docs
ev
class_one
Psi
virtual_zero_section
derivated_functor
reduced_virtual_zero_section
```