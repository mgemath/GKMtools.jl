# Generalities on GKM graphs

## Definition

GKM spaces have been introduced in [GKM98](@cite). Many of the combinatorial definitions in this package follow [GZ98](@cite). For the purpose of this package, by *GKM variety* we mean a smooth projective varieties over $\mathbb{C}$ with an algebraic torus action such that the action has a finite number of fixed points and a finite number of 1-dimensional orbits.

The *GKM graph* associated to a torus $T$ acting on a GKM variety $X$ is the following datum:
* A graph having the fixed points as vertices, such that two vertices are connected by an unoriented) edge if there is a 1-dimensional orbit passing through the two fixed points.
* An axial function $\mathrm{w}\colon E \rightarrow M$ from the set of oriented edges of the graph to the weight lattice $M$ of $T$. (By *oriented edge* we mean an unoriented edge of the graph plus a choice of orientation.)

In this package, the codomain $M$ of the axial function can be a free $\mathbb{Z}$-module or a free $\mathbb{Q}$-module.
Since the GKM graph of a GKM variety is always regular (with the valency of every vertex being the complex dimension of the space), this package assumes that GKM graphs are regular.

Famous examples of GKM varities include projective space, (generalised/partial) flag varieties, smooth Schubert varieties, and smooth toric varieties, see [Standard Constructions](STDconstructions.md).

This package represents GKM varieties purely by their GKM graphs. For some applications, the additional datum of a *GKM connection* is necessary, see [Connections](Connections.md).

## Index

```@index
```