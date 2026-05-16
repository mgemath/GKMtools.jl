# Cohomology

Let $X$ be a GKM space with respect to the complex torus $T$.
By [GKM98](@cite) and [GZ01; Theorem 1.7.3](@cite), we have the following description of its equivariant cohomology ring $H^*_T(X;\mathbb{Q})$.

Each fixed point $x\in X$ gives a ring homomorphism $H^*_T(X;\mathbb{Q})\rightarrow H^*_T(\{x\};\mathbb{Q})\cong \mathbb{Q}[\mathfrak{t}]$.
We may combine these maps by taking all fixed points at once.
Let $V$ be the set of fixed points of $X$ and $E$ the set of $T$-invariant rational curves.
That is, $V$ are the vertices of the GKM graph and $E$ the edges.
Then the map

```math
    H^*_T(X;\mathbb{Q}) \longrightarrow \bigoplus_{x\in V} \mathbb{Q}[\mathfrak{t}]
```

is injective and its image consists of all $(f_x)_{x\in V}$ such that $f_{\text{src}(e)}\equiv f_{\text{dst}(e)}$ mod $w(e)$ for all edges $e\in E$, where $w(e)$ is the weight of $e$.

We may further identify $\mathbb{Q}[\mathfrak{t}]\cong \mathbb{Q}[t_1,\dots,t_r]$ where $r=\dim_{\mathbb{C}}(T)$.
Hence, this package represents elements of $H^*_T(X;\mathbb{Q})$ as tuples of polynomials indexed by vertices of the GKM graph.

```@docs
is_gkm_class
point_class
poincare_dual
weight_class
euler_class
integrate_gkm_class
integrate
first_chern_class
chern_class
```