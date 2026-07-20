# Orbifold GKM graphs

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

## Inertia stack

```@docs
inertia_stack
twisted_sectors
untwisted_sector
sector_count
```

## Orbifold vector bundles

```@docs
orbifold_vector_bundle
orbifold_line_bundle
fiber_representation
```
