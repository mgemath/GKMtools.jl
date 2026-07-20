# Generalities on GKM graphs

## Definition

!!! note
    Although the definition below is written for smooth algebraic GKM spaces, the functionality of this package applies equally well to Hamiltonian GKM spaces as discussed in the supporting article of this package.
    By *GKM space* we will mean either a GKM variety as defined below, or a Hamiltonian GKM space.

GKM spaces have been introduced in [GKM98](@cite). Many of the combinatorial definitions in this package follow [GZ98](@cite). For the purpose of this package, by *GKM variety* we mean a smooth projective variety over ``\mathbb{C}`` with an algebraic torus action such that the action has a finite number of fixed points and a finite number of 1-dimensional orbits.

The *GKM graph* associated to a torus ``T`` acting on a GKM variety ``X`` is the following datum:

* A graph having the fixed points as vertices, such that two vertices are connected by an unoriented edge if there is a 1-dimensional orbit passing through the two fixed points.
* An axial function ``\mathrm{w}\colon E \rightarrow M`` from the set of oriented edges of the graph to the weight lattice ``M`` of ``T``. By *oriented edge* we mean an unoriented edge of the graph plus a choice of orientation.

In this package, the codomain ``M`` of the axial function can be a free ``\mathbb{Z}``-module or a free ``\mathbb{Q}``-module.
Since the GKM graph of a GKM variety is always regular, with the valency of every vertex being the complex dimension of the space, this package assumes that GKM graphs are regular.

Famous examples of GKM varieties include projective space, generalized/partial flag varieties, smooth Schubert varieties, and smooth toric varieties; see [Standard Constructions](STDconstructions.md).

This package represents GKM varieties purely by their GKM graphs. For some applications, the additional datum of a *GKM connection* is necessary; see [Connections](Connections.md).

!!! note
    We have added support for non-compact GKM spaces, which arise for example from quasi-projective algebraic GKM spaces or as total spaces of GKM vector bundles over GKM spaces. On the level of GKM graphs, this means that *standalone flags* (sometimes called *semi-infinite edges*) are allowed:

    * Each vertex ``p`` of the GKM graph has a set of flags. These are given by the ``T``-invariant linear subspaces of ``T_pX``.
    * The axial function assigns to each flag the ``T``-weight of that linear subspace.
    * Two flags at different vertices form an edge if and only if they correspond to tangent spaces of the same 1-dimensional orbit.
    * Every edge consists of precisely two flags.

## Orbifold GKM graphs

Orbifold GKM graphs encode the same fixed-point and 1-dimensional-orbit combinatorics, together with the finite stabilizer data that appears in orbifold charts.
In this package an orbifold GKM graph is represented by:

* An underlying GKM graph with vertices, flags, edges, and axial weights.
* For every vertex ``p``, a finite abelian isotropy group ``G_p`` and its representation on the tangent space ``T_pX``.
* For every flag, a finite abelian isotropy group of the corresponding 1-dimensional stratum.
* For every flag incident to a vertex, an embedding from the flag isotropy group into the vertex isotropy group.

The isotropy groups are stored by their cyclic factors, for example ``[2, 3]`` represents ``\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/3\mathbb{Z}``.
The tangent representation at a vertex is stored as an integer matrix whose rows correspond to cyclic factors of the vertex isotropy group and whose columns correspond to tangent weight lines.
Similarly, the flag embedding matrix records how the stabilizer of a 1-dimensional stratum maps into the stabilizer of an endpoint.

Geometrically, this extra data allows the package to distinguish an ordinary smooth GKM graph from a stacky or orbifold one, even when the underlying graph and axial weights look similar.
For instance, weighted projective spaces have orbifold GKM graphs: the vertices are the torus fixed points, while some vertices and 1-dimensional strata may carry non-trivial finite stabilizers.

This orbifold data is used by constructions such as the inertia stack, orbifold vector bundles, and Chern classes of orbifold GKM vector bundles.

## Index

```@index
```
