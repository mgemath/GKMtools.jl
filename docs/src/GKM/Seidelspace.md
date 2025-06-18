# Seidel Space

We follow [Iri17; Definition 3.2](@cite) in our description of the Seidel space.
Let $X$ be a (smooth) GKM variety with torus action by $T$.
Let $\iota\colon\mathbb{C}^\times\rightarrow T$ be a group homomorphism.

The *Seidel space* $S_X$ associated to this datum is a $X$-bundle over $\mathbb{P}^1$ which is trivialized over
```math
  \begin{align*}
    U_0 &= \{[1:u] : u\in\mathbb{C}\}\subset\mathbb{P}^1\\
    U_\infty &= \{[v:1] : v\in\mathbb{C}\}\subset\mathbb{P}^1.
  \end{align*}
```
The transition function on $U_0\cap U_\infty$ is given by

```math
  U_0\times X \ni ([1:u], x) \longmapsto ([u^{-1}:1], \iota(u)\cdot x) \in U_\infty\times X.
```

Globally, it can be described as
```math
  S_X = \left( \mathbb{C}^2\setminus\{0\}\times X \right) / \mathbb{C}^\times,
```
where $\lambda\in\mathbb{C}^\times$ acts by $\lambda\cdot ((v, u), x) := ((\lambda v, \lambda u), \iota(\lambda)\cdot x)$.


Let $\widehat{T}:= T\times\mathbb{C}^\times$.
Then $\widehat{T}$ acts on $S_X$ by $(t, s)\cdot ([(v, su), t\cdot x])$.

That is, $\widehat{T}$ acts on $S_X\rightarrow\mathbb{P}^1$ via the $T$-action on the fibres $X$ and via the standard $\mathbb{C}^\times$-action on the base $\mathbb{P}^1$.

!!! note
    - If $X$ is a GKM space with respect to $T$, then $S_X$ is a GKM space with respect to $\widehat{T}$, which is implemented in this package.
    - The GKM graph and connection of $S_X$ is described in the supporting paper.

```@docs
Seidel_space
```