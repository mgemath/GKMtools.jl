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

If $G$ is the GKM graph of a GKM variety $X$, then these integers $a$ are the degrees of the equivariant line bundles into which $TX$ splits when restricted to the invariant rational curve represented by $e$.

## Existence and uniqueness of connections

Given a GKM graph $G$ that comes from a GKM variety $X$, it always has a connection for the geometric reason sketched above. However, ``G`` may admit more than one connection.

The following are sufficient conditions for the existence of a unique connection of ``G``.
 * The valency of $G$ is at least 3 and $G$ is $3$-independent, i.e. the weights of every three edges starting at the same vertex are linearly independent.
 * The valency of $G$ is at most 2.

In those cases, the connection can be computed using `get_connection`.

If neither of these two conditions hold and $G$ is not the output of a standard construction, a choice of connection can be specified manually using `set_connection!`.

```@docs
get_connection(::GKMtools.GKM_graph)
get_any_connection(::GKMtools.GKM_graph)
build_gkm_connection
set_connection!
isvalid_connection
is_compatible_with_connection
```