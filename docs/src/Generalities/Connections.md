# Connections

Let ``e`` be a directed edge ``p\rightarrow q``, and let ``F_p`` (resp., ``F_q``) be the set of all flags at ``p`` (resp., ``q``). Following [GZ98](@cite), a *connection* along ``e`` is a bijection

```math
\nabla_e\colon F_p \longrightarrow F_q.
```

A connection ``\nabla`` of a GKM graph ``G`` is a family of connections ``\nabla = \{\nabla_e \}_{e \in E}``, where ``E`` is the set of all edges of ``G``, such that ``\nabla_{-e}=\nabla_{e}^{-1}``.

The connection is compatible with the axial function ``\mathrm{w}`` of ``G`` if for all ``e' \in F_p``, there exists a coefficient ``a`` depending on ``e`` and ``e'`` such that

```math
\mathrm{w}(\nabla_e(e')) = \mathrm{w}(e') - a \mathrm{w}(e).
```

If $G$ is the GKM graph of a GKM space $X$, then these integers $a$ are the degrees of the equivariant line bundles into which $TX$ splits when restricted to the invariant rational curve represented by $e$.

## Existence and uniqueness of connections

Given a GKM graph $G$ that comes from a GKM space $X$, it always has a connection for the geometric reason sketched above. However, ``G`` may admit more than one connection.

The following are sufficient conditions for the existence of a unique connection of ``G``.
 * The valency of $G$ is at least 3 and $G$ is $3$-independent, i.e. the weights of every three edges starting at the same vertex are linearly independent.
 * The valency of $G$ is at most 2.

Toric varieties always satisfy those conditions, hence they admits a unique connection.

Some other varieties admits more than connection, for example 'generalized_gkm_graph'.

When a GKM graph is created, unless otherwise set, a connection is always created using a combinatorial algorithm. This connection is created using the following function:

```@docs
build_gkm_connection
```

In order to show the connection of a GKM graph, use the following function:

```@docs
print_connection
```

Some GKM varieties allow the choice of more GKM connections, for example the variety of complete flags in $\mathbb{C}^3$, see [Standard Constructions](STDconstructions.md).

```jldoctest conn
julia> R = root_system(:A, 2)
Root system of rank 2
  of type A2

julia> G1 = generalized_gkm_flag(R)
GKM graph with 6 nodes, valency 3 and axial function:
s1 -> id => (-1, 1, 0)
s2*s1 -> s1 => (0, -1, 1)
s1*s2*s1 -> id => (-1, 0, 1)
s1*s2*s1 -> s2*s1 => (-1, 1, 0)
s2 -> id => (0, -1, 1)
s2 -> s2*s1 => (1, 0, -1)
s1*s2 -> s1 => (-1, 0, 1)
s1*s2 -> s1*s2*s1 => (0, 1, -1)
s1*s2 -> s2 => (-1, 1, 0)
Birkhoff-Grothendieck connection for GKM graph with 6 nodes and valency 3
```
By default, `generalized_gkm_flag` equips the graph with the connection coming from the Birkhoff-Grothendieck decomposition. In order to have the connection using the root data, use the following:

```jldoctest conn
julia> G2 = generalized_gkm_flag(R; connection = :cartan)
GKM graph with 6 nodes, valency 3 and axial function:
s1 -> id => (-1, 1, 0)
s2*s1 -> s1 => (0, -1, 1)
s1*s2*s1 -> id => (-1, 0, 1)
s1*s2*s1 -> s2*s1 => (-1, 1, 0)
s2 -> id => (0, -1, 1)
s2 -> s2*s1 => (1, 0, -1)
s1*s2 -> s1 => (-1, 0, 1)
s1*s2 -> s1*s2*s1 => (0, 1, -1)
s1*s2 -> s2 => (-1, 1, 0)
Cartan connection for GKM graph with 6 nodes and valency 3
```
One can use also the algorithm connection:

```jldoctest conn
julia> G3 = generalized_gkm_flag(R; connection = :algorithm);
```
In this example, the Birkhoff-Grothendieck connection and the Cartan connection do not coincide.
```jldoctest conn
julia> print_connection(G1)
Birkhoff-Grothendieck connection for GKM graph with 6 nodes and valency 3
Connection:
s2 -> id => [2, 3, 1]
s1*s2*s1 -> s1*s2 => [1, 3, 2]
s2*s1 -> s1*s2*s1 => [3, 1, 2]
s1*s2*s1 -> id => [3, 2, 1]
s2 -> s1*s2 => [1, 2, 3]
s1 -> id => [1, 3, 2]
id -> s2 => [3, 1, 2]
id -> s1 => [1, 3, 2]
s2*s1 -> s1 => [3, 1, 2]
s2*s1 -> s2 => [1, 2, 3]
s2 -> s2*s1 => [1, 2, 3]
id -> s1*s2*s1 => [3, 2, 1]
s1 -> s1*s2 => [3, 1, 2]
s1 -> s2*s1 => [2, 3, 1]
s1*s2*s1 -> s2*s1 => [2, 3, 1]
s1*s2 -> s1*s2*s1 => [1, 3, 2]
s1*s2 -> s2 => [1, 2, 3]
s1*s2 -> s1 => [2, 3, 1]
a_i's:
s2 -> id => [2, 0, 0]
s1*s2*s1 -> s1*s2 => [0, 0, 2]
s2*s1 -> s1*s2*s1 => [0, 0, 2]
s1*s2*s1 -> id => [2, 1, 1]
s2 -> s1*s2 => [1, 1, 2]
s1 -> id => [2, 0, 0]
id -> s2 => [0, 2, 0]
id -> s1 => [2, 0, 0]
s2*s1 -> s1 => [2, 1, 1]
s2*s1 -> s2 => [0, 2, 0]
s2 -> s2*s1 => [0, 2, 0]
id -> s1*s2*s1 => [1, 1, 2]
s1 -> s1*s2 => [0, 2, 0]
s1 -> s2*s1 => [1, 1, 2]
s1*s2*s1 -> s2*s1 => [0, 2, 0]
s1*s2 -> s1*s2*s1 => [0, 2, 0]
s1*s2 -> s2 => [1, 1, 2]
s1*s2 -> s1 => [2, 0, 0]

julia> print_connection(G2)
Cartan connection for GKM graph with 6 nodes and valency 3
Connection:
s2 -> id => [2, 1, 3]
s1*s2*s1 -> s1*s2 => [3, 1, 2]
s2*s1 -> s1*s2*s1 => [1, 3, 2]
s1*s2*s1 -> id => [3, 2, 1]
s2 -> s1*s2 => [1, 2, 3]
s1 -> id => [1, 2, 3]
id -> s2 => [2, 1, 3]
id -> s1 => [1, 2, 3]
s2*s1 -> s1 => [3, 1, 2]
s2*s1 -> s2 => [3, 2, 1]
s2 -> s2*s1 => [3, 2, 1]
id -> s1*s2*s1 => [3, 2, 1]
s1 -> s1*s2 => [2, 1, 3]
s1 -> s2*s1 => [2, 3, 1]
s1*s2*s1 -> s2*s1 => [1, 3, 2]
s1*s2 -> s1*s2*s1 => [2, 3, 1]
s1*s2 -> s2 => [1, 2, 3]
s1*s2 -> s1 => [2, 1, 3]
a_i's:
s2 -> id => [2, -1, 1]
s1*s2*s1 -> s1*s2 => [1, -1, 2]
s2*s1 -> s1*s2*s1 => [1, -1, 2]
s1*s2*s1 -> id => [2, 1, 1]
s2 -> s1*s2 => [1, 1, 2]
s1 -> id => [2, -1, 1]
id -> s2 => [-1, 2, 1]
id -> s1 => [2, -1, 1]
s2*s1 -> s1 => [2, 1, 1]
s2*s1 -> s2 => [1, 2, -1]
s2 -> s2*s1 => [-1, 2, 1]
id -> s1*s2*s1 => [1, 1, 2]
s1 -> s1*s2 => [-1, 2, 1]
s1 -> s2*s1 => [1, 1, 2]
s1*s2*s1 -> s2*s1 => [1, 2, -1]
s1*s2 -> s1*s2*s1 => [-1, 2, 1]
s1*s2 -> s2 => [1, 1, 2]
s1*s2 -> s1 => [2, -1, 1]

julia> print_connection(G3, verbose = false)
Algorithmic connection for GKM graph with 6 nodes and valency 3
```
