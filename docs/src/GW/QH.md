# Quantum Cohomology

Much of this can be found in [CK99; Chapter 8 and 9.3](@cite).
Let $X$ be a GKM space (see [Definition](../GKM/GKM.md#Definition)).
Its (small) *equivariant quantum cohomology* $QH_T^*(X)$ is given additively by
$H_T^*(X;\mathbb{Q})\otimes \widehat{\mathbb{Q}[H_2^\text{eff}(X;\mathbb{Z})]}$, where $\widehat{\mathbb{Q}[H_2^\text{eff}(X;\mathbb{Z})]}$ is the completion
of the semigroup ring $H_2^\text{eff}(X;\mathbb{Z})$ of effective curve classes.
The element corresponding to $\beta\in H_2^\text{eff}(X;\mathbb{Z})$ is written as $q^\beta$.

The $H_T(\text{pt};\mathbb{Q})$-module $QH_T^*(X)$ is a commutative associative unital $H_T(\text{pt};\mathbb{Q})$-algebra via the (small) *equivariant quantum product* $\ast$ defined as follows.
For every classes $a,b,c\in H_T^*(X;\mathbb{Q})$ we have
```math
  \langle a \ast b, c \rangle = \sum_{\beta\in H_2^\text{eff}(X;\mathbb{Z})} GW^T_{X,0,3,\beta}(a,b,c) \cdot q^\beta
```
where:
 * The equivariant Poincaré pairing is given by $\langle a,b\rangle := \int_X a\cup b\in H_T^*(\text{pt};\mathbb{Q})$, where we use equivariant integration,
 * We denote by $GW^T_{X,0,3,\beta}(a,b,c)\in H_T^*(\text{pt};\mathbb{Q})$ the equivariant Gromov--Witten invariant for $X$ in class $\beta$ of genus $0$ with $3$ marked points.

Note that setting all the equivariant parameters $t_1,\dots,t_{\dim_\mathbb{C}(T)}$ to zero recovers the standard (small, non-equivariant) quantum product.

## Structure Constants

Given any linear basis $e_1,\dots,e_N$ of $H_T^*(X;\mathbb{Q})\otimes \text{Frac}(H_T^*(\text{pt}; \mathbb{Q}))$ over $\text{Frac}(H_T^*(\text{pt}; \mathbb{Q}))\cong \mathbb{Q}(t_1,\dots,t_r)$, the associated *structure constants* are
$(c_{i,j,k})_{i,j,k=1}^N$ where
```math
c_{i,j,k}\in \mathbb{Q}(t_1,\dots,t_r)\otimes \widehat{\mathbb{Q}[H_2^\text{eff}(X;\mathbb{Z})]}
```
is the coefficient of $e_k$ in $e_i\ast e_j$.
Let us denote the coefficient of $q^\beta$ in $c_{i,j,k}$ by $c_{i,j,k}(\beta)\in\mathbb{Q}(t_1,\dots,t_r)$.

!!! note
    When $(X,T)$ is a compact Hamiltonian or projective GKM space, and $e_1,\dots,e_N$ is chosen to be a linear basis of $H_T^*(X;\mathbb{Q})$ over $H_T^*(\text{point};\mathbb{Q})\cong \mathbb{Q}[t_1,\dots,t_r]$, we have
    ```math
        c_{i,j,k}(\beta)\in\mathbb{Q}[t_1,\dots,t_r]
    ```
    for all $i,j,k,\beta$.

### Choice of basis

Since structure constants depend on a choice of basis, let us introduce some common choices and explain structure constants for [mixed bases](#Mixing-the-bases).

#### The standard basis

In general, when no basis is specified explicitly, all structure constants are computed with respect to the basis $e_1,\dots,e_N$, which we call the *standard basis*.
Namely, let $N$ be the number of vertices of the GKM graph $G$. For each $i\in\{1,\dots,N\}$, let
```math
    e_i \in H_T^*(X;\mathbb{Q})\otimes \mathbb{Q}(t_1,\dots,t_r)
```
be class that localizes to $1$ at the $i$-th fixed point and to $0$ at every other fixed point.
This uniquely defines a class by the GKM theorem [GKM98](@cite).

!!! note
    - The $e_i$ are only well-defined as elements of $H_T^*(X;\mathbb{Q})\otimes \mathbb{Q}(t_1,\dots,t_r)$, not of $H_T^*(X;\mathbb{Q})$.
    - Since the standard basis over $\mathbb{Q}(t_1,\dots,t_r)$, the resulting structure constants $c_{i,j,k}(\beta)$ can fail to be polynomials in $t_1,\dots,t_r$ even when $G$ is the GKM graph of a compact Hamiltonian or projective GKM space $(X,T)$.

In our implementation, $e_i$ is printed as `e[i]` (see [Cohomology](../GKM/Cohomology.md)).

#### The fixed point basis

The *fixed point basis* is given by $f_1,\dots,f_N$, where $N$ is the number of vertices of $G$ and
```math
    f_i = \left( \prod_{\epsilon\in E(G)_i} \alpha(\epsilon) \right) e_i \in H_T^*(X;\mathbb{Q}).
```
Mathematically, we have $\prod_{\epsilon\in E(G)_i} \alpha(\epsilon) = e_T(T_{p_i}X)$ when $G$ is realized by the GKM space $X$ and $p_i\in X^T$ corresponds to the $i$-th vertex of $G$.

!!! note
    - Sometimes, we use the notation $f_i = PD(p_i)$ as $f_i$ is the *equivariant Poincaré dual* of $p_i\in X^T$.
    - The collection $(f_i)_{i=1}^N$ is not a $\mathbb{Q}[t_1,\dots,t_r]$-linear basis for $H_T^*(X;\mathbb{Q})$, but it is a $\mathbb{Q}(t_1,\dots,t_r)$-linear basis for $H_T^*(X;\mathbb{Q})\otimes \mathbb{Q}(t_1,\dots,t_r)$.
    - The element $f_i$ can be obtained computationally as `point_class(i, G)` (see [`point_class`](@ref)).

#### Mixing the bases
In some functions (such as [`QH_structure_constants`](@ref)) we use two different bases to define $c_{i,j,k}(\beta)$.
Namely, given two bases $(a_i)$ and $(b_i)$, we let $c_{i,j,k}(\beta)$ be the coefficient of $b_k q^\beta$ in $a_i \ast a_j$.
In this case, we say that $(c_{i,j,k}(\beta))$ are the structure constants with respect to the *input basis* $(a_i)$ and *output basis* $(b_i)$.

```@docs
QH_structure_constants
QH_structure_constants_in_basis
QH_supporting_curve_classes
```

## Quantum Arithmetic

The implementation represents an element as a finite dictionary
``\beta \mapsto a_\beta`` corresponding to
``\sum_\beta a_\beta q^\beta``. Coefficients are localized GKM classes.
Because a general GKM space can have infinitely many contributing effective
classes, a parent records an explicit finite set of Novikov degrees. Products
are projected to this truncation.

```jldoctest quantum_manual
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta = GKMtools.curve_class(P1, Oscar.Edge(1, 2));

julia> QH = small_equivariant_quantum_cohomology(P1; degrees=[beta]);

julia> p = quantum_class(QH, point_class(P1, 1));

julia> q = quantum_class(QH, 1; degree=beta);

julia> iszero(coefficient(p, beta)) && isone(coefficient(q, beta))
true

julia> p * one(QH) == p
true
```

Degree-zero multiplication is the ordinary equivariant cup product. For a
positive effective class, [`quantum_product`](@ref) computes the coefficient
with [`gromov_witten_nomarks`](@ref), integrating the factors and fixed-point
classes over the image curve. Multiplication of quantum classes combines these
coefficients and retains only the degrees in the parent truncation.

## API

```@docs
SmallEquivariantQuantumCohomology
SmallEquivariantQuantumClass
small_equivariant_quantum_cohomology
truncation_degrees
quantum_class
QH_class
quantum_coefficients
coefficient
quantum_product
small_quantum_cohomology_ring
```

The standard operations `+`, `-`, `*`, `^`, `zero`, and `one` are supported.

The constructor [`small_quantum_cohomology_ring`](@ref) packages a finite set
of computed quantum corrections into an Oscar graded quotient. It returns the
quotient, its cohomology basis grouped by codimension, and its Novikov
generators.

```@setup
#= 
<!-- ```@docs
QH_class
*(::GKMtools.QHRingElem, ::GKMtools.QHRingElem)
quantum_product
quantum_product_at_q1
```
=#
#=
Quantum product with $c_1^T(TX)$
```@docs
c1_at_q1
conjecture_O_eigenvalues
```
=#
#=
Twisted versions
```@docs
twisted_c1_matrix
twisted_c1_matrix_at_q1
```
=#
#=
Sanity checks

```@docs
QH_is_commutative
QH_is_associative
QH_is_homogeneous
QH_is_polynomial
``` -->
=#
```
