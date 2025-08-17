# Quantum Cohomology

Much of this can be found in [CK99; Chapter 8 and 9.3](@cite).
Let $X$ be a (smooth projective) GKM variety.
Its (small) *equivariant quantum cohomology* $QH_T^*(X)$ is given additively by
$H_T^*(X;\mathbb{Q})\otimes \widehat{\mathbb{Q}[H_2^\text{eff}(X;\mathbb{Z})]}$, where $\widehat{\mathbb{Q}[H_2^\text{eff}(X;\mathbb{Z})]}$ is the completion
of the semigroup ring $H_2^\text{eff}(X;\mathbb{Z})$ of effective curve classes.
The element corresponding to $\beta\in H_2^\text{eff}(X;\mathbb{Z})$ is written as $q^\beta$.

The $H_T(\text{pt};\mathbb{Q})$-module $QH_T^*(X)$ is a commutative associative unital $H_T(\text{pt};\mathbb{Q})$-algebra via the (small) *equivariant quantum product* $\ast$ defined as follows.
For every classes $a,b,c\in H_T^*(X;\mathbb{Q})$ we have
```math
  \langle a \ast b, c \rangle = \sum_{\beta\in H_2^\text{eff}(X;\mathbb{Z})} GW^T_{0,3,\beta}(a,b,c) \cdot q^\beta
```
where:
 * The equivariant Poincaré pairing is given by $\langle a,b\rangle := \int_X a\cup b\in H_T^*(\text{pt};\mathbb{Q})$, where we use equivariant integration,
 * We denote by $GW^T_{0,3,\beta}(a,b,c)\in H_T^*(\text{pt};\mathbb{Q})$ the equivariant Gromov--Witten invariant for $X$ in class $\beta$ of genus $0$ with $3$ marked points.

Note that setting all the equivariant parameters $t_1,\dots,t_{\dim_\mathbb{C}(T)}$ to zero recovers the standard (small, non-equivariant) quantum product.

