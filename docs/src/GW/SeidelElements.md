# Seidel elements / Shift operators

This section deals with certain equivariant Gromov--Witten invariants on [Seidel spaces](../GKM/Seidelspace.md).
The definition used in this package should be carefully compared to [Iri17; Section 3](@cite) and [MO12; Chapter 8](@cite) *(shift operators)*.
A (non-equivariant) symplectic account can be found in [MS12; Section 11.4](@cite) *(Seidel representation)*.

Let $X$ be a (smooth projective) GKM variety with torus action by $T$, and let $\iota\colon  \mathbb{C}^\times \rightarrow T$ be a group homomorphism.
This gives rise to the $X$-bundle $\pi\colon S_X\rightarrow \mathbb{P}^1$, where $S_X$ is the [Seidel space](../GKM/Seidelspace.md) associated to $(X, \iota)$.
Recall that $S_X$ is a GKM space with respect to $\widehat{T}:= T\times\mathbb{C}^\times$, where the extra copy of $\mathbb{C}^\times$ acts by rotating the
base $\mathbb{P}^1$.

Let $(e_i)$ and $(e^i)$ be dual bases of $H_T^*(X;\mathbb{Q})$ with respect to the $T$-equivariant Poincaré pairing on $X$.
Note also that we have a $H_T^*(\text{pt};\mathbb{Q})$-linear pushforward map $(i_\infty)_*\colon H_T^*(X;\mathbb{Q}) \rightarrow H_{\widehat{T}}^*(S_X;\mathbb{Q})$ raising degree by one.

Finally, let $H_2^\text{sec}(S_X;\mathbb{Z})$ be the set of effective section curve classes in $S_X$, i.e.,
effective curve classes that project to $[\mathbb{P}^1]$ under $\pi\colon S_X\rightarrow\mathbb{P}^1$.
The additive group $H_2^\text{sec}(S_X;\mathbb{Z})$ is an $H_2^\text{eff}(X;\mathbb{Z})$-torsor, so after picking some $\beta_0\in H_2^\text{sec}(S_X;\mathbb{Z})$ there is an identification
$r\colon H_2^\text{sec}(S_X;\mathbb{Z})\stackrel{\cong}{\longrightarrow} H_2^\text{eff}(X;\mathbb{Z})$ that sends $\beta_0$ to $0$.

The *(equivariant) Seidel element* associated to $\iota$ is

```math
  \mathcal{S}(\iota) := \sum_{\beta\in H_2^\text{sec}(S_X;\mathbb{Z})} \left[GW^{\widehat{T}}_{0,1,\beta}((i_\infty)_*(e_i))\right]_{\hat{t} = 0} e^i q^{r(\beta)} \in QH_T^*(X)
```
where $\hat{t}$ is the equivariant parameter for the extra $\mathbb{C}^\times$ in $\widehat{T}$.

!!! note
    Let $\Lambda:= \text{Hom}(\mathbb{C}^\times, T)$ be the cocharacter lattice of $T$.
    It is a key property of the equivariant Seidel elements that the map
    ```math
      \mathcal{S}\colon \Lambda \rightarrow QH_T^*(X)^\times
    ```
     is a group homomorphism, where $QH_T^*(X)^\times$ is endowed with the equivariant quantum product.

```@docs
Seidel_element
```