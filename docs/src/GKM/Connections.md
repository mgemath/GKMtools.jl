# Connections

Let ``e`` be an directed edge ``p\rightarrow q``, and let ``E_p`` (resp., ``E_q``) be the set of all directed edges starting from ``p`` (resp., ``q``). Following [GZ98](@cite), a *connection* along ``e`` is a bijection

```math
\nabla_e\colon E_p \longrightarrow E_q.
```

A connection ``\nabla`` of a GKM graph ``G`` is a family of connections ``\nabla = \{\nabla_e \}_{e \in E}``, where ``E`` is the set of all edges of ``G``, such that ``\nabla_{-e}=\nabla_{e}^{-1}``.

The connection is compatible with the axial function ``\mathrm{w}`` of ``G`` if for all ``e' \in E_p``, there exists an integer ``a`` depending on ``e`` and ``e'`` such that

```math
\mathrm{w}(\nabla_e(e')) = \mathrm{w}(e') - a \mathrm{w}(e).
```

If $G$ is the GKM graph of a GKM space $X$, then these integers $a$ are the degrees of the equivariant line bundles into which $TX$ splits when restricted to the invariant rational curve represented by $e$.

## Existence and uniqueness of connections

Given a GKM graph $G$ that comes from a GKM space $X$, it always has a connection for the geometric reason sketched above. However, ``G`` may admit more than one connection.

The following are sufficient conditions for the existence of a unique connection of ``G``.
 * The valency of $G$ is at least 3 and $G$ is $3$-independent, i.e. the weights of every three edges starting at the same vertex are linearly independent.
 * The valency of $G$ is at most 2.

In those cases, the connection can be computed using `get_connection`.

If neither of these two conditions hold and $G$ is not the output of a standard construction, a choice of connection can be specified manually using `set_connection!`.

## Support for standalone flags

!!! note
    In our ongoing efforts to support GKM graphs with standalone flags, a connection is no longer represented
    as a dictionary from pairts of edges to edges, but as a bijection of flag indices for each edge.
    That is, for every oriented edge $e\in E(G)^\pm$, we represent $\nabla_e$ as `Vector{Int64}`.
    If entry `i` of this vector is `j` then $\nabla_e$ sends flag `i` at `src(e)` to flag `j` at `dst(e)`.

```@docs
get_connection(::GKMtools.AbstractGKM_graph)
get_any_connection(::GKMtools.AbstractGKM_graph)
build_GKM_connection
set_connection!
isvalid(::GKMtools.GKM_connection; ::Bool)
is_compatible_with_connection
```