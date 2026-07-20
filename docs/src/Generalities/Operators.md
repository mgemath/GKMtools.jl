# Operators

## GKM Subgraphs
The following figure shows the GKM subgraph of a Schubert variety $X_w$ in the full flag variety of $\mathbb{C}^3$.
Edge labels illustrate the axial function, while the vertex labels describe the equivariant Poincaré dual of $X_w$ via localization.
![Illustration of GKM subspaces](../img/subgraph.svg)
```@docs
subgraph_from_vertices
ambient_graph
subgraph
vertex_to_ambient
ambient_to_vertex
poincare_dual
```

## Blowups
The following figure illustrates the effect of blowups along a GKM subgraph (red) on the underlying graph.
![Illustration of blowups along sub-GKM-graphs](../img/blowup.svg)
```@docs
blow_up
weighted_blow_up
```