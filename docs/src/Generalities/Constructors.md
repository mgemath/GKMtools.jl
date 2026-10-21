```@meta
DocTestSetup = :(using Oscar, GKMtools)
```

# Constructors

The main way to construct GKM graphs is using one of the constructors provided in the page [Examples](../Examples/Examples.md). The construction of a GKM graph from scratch is possible but not advised.

# Example
Here we construct the $2$-dimensional projective space from scratch. 
```jldoctest scratch
julia> using GKMtools, Oscar
```
Let us consider the coordinate torus $T=(\mathbb{C}^*)^3$ acting on $\mathbb{P}^2$:
```math
(t_0, t_1, t_2) \cdot [x_0:x_1:x_2] = [t_0x_0:t_1x_1:t_2x_2]
```

The fixed points are
```math
p_0 = [1:0:0], p_1 = [0:1:0], p_2 = [0:0:1].
```

Let us construct the character lattice of $T$. It is a Oscar free module. 
```jldoctest scratch
julia> M = free_module(ZZ, 3) # also possible to use QQ, but we want to see the integral structure
Free module of rank 3 over ZZ

julia> t_0, t_1, t_2 = gens(M)
3-element Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}:
 (1, 0, 0)
 (0, 1, 0)
 (0, 0, 1)

julia> characters = (t_0, t_1, t_2)
((1, 0, 0), (0, 1, 0), (0, 0, 1))
```
We have $3$ fixed points and the one-dimensional $T$-orbits connect every pair of fixed points,
so the underlying GKM graph is a triangle.
```jldoctest scratch
julia> g = Graph{Undirected}(3); # Oscar 3-vertices graph with no edges

julia> add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 2, 3);
```
Let us construct the vertices of our GKM graph. Our vertex contains a string which is its label. The label appears when we print the GKM graph, it does not have any other goal. Two vertices cannot have the same label.

```jldoctest scratch
julia> vertex_labels = Vertex.([
           "p_0 = [1:0:0]",
           "p_1 = [0:1:0]",
           "p_2 = [0:0:1]",
       ])
3-element Vector{Vertex}:
 Vertex("p_0 = [1:0:0]")
 Vertex("p_1 = [0:1:0]")
 Vertex("p_2 = [0:0:1]")
```
Let us now define the values of the axial map. For each vertex `p`, we define a vector `A` of `FlagWeight{ZZRingElem}`, which will contain the value of the axial map. In such a way that `A[i]` will give the value of the `i`-flag starting from the vertex `p`.  That will be nothing more that an element of $M$. We start by constructing a vector of length equal to the number of vertices, storing all those vectors `A`.
```jldoctest scratch
julia> flags = [
           FlagWeight{ZZRingElem}[]
           for _ in 1:3
       ]
3-element Vector{Vector{FlagWeight{ZZRingElem}}}:
 []
 []
 []
```
At any point $p_i$, the tangent weight pointing toward $p_j$ is $t_j-t_i$. 
For each edge `e`, `edge_flags[e]=(a, b)` means that this edge defines the flag number `a` of the vertex `src(e)` and the flag number `b` of the vertex `dst(e)`.
```jldoctest scratch
julia> edge_flags = Dict{Edge,Tuple{Int,Int}}();

julia> for e in edges(g)
           i = src(e)
           j = dst(e)
           A_i = characters[j] - characters[i]
           # Opposite orientations of the same invariant P¹ have
           # opposite tangent weights.
           push!(flags[i], FlagWeight{ZZRingElem}(A_i))
           push!(flags[j], FlagWeight{ZZRingElem}(-A_i))
           edge_flags[e] = (
               length(flags[i]),
               length(flags[j]),
           )
       end
```
Assemble the raw combinatorial GKM data.
```jldoctest scratch
julia> core_data = GKMCombinatorialData{
           ZZRingElem,
           Vertex,
           FlagWeight{ZZRingElem},
       }(
           g,
           M,
           vertex_labels,
           flags,
           edge_flags,
       )
GKM Combinatorial Data
  Graph: Undirected graph with 3 nodes and 3 edges
  Character Lattice: Free module of rank 3 over ZZ
  Vertex Labels: ["p_0 = [1:0:0]", "p_1 = [0:1:0]", "p_2 = [0:0:1]"]
  Flag weights defined at 3 vertices

```
Finally, construct the GKM graph. A connection will be given automatically. 
```jldoctest scratch
julia> P2 = gkm_graph(core_data)
GKM graph with 3 nodes, valency 2 and axial function:
p_1 = [0:1:0] -> p_0 = [1:0:0] => (1, -1, 0)
p_2 = [0:0:1] -> p_0 = [1:0:0] => (1, 0, -1)
p_2 = [0:0:1] -> p_1 = [0:1:0] => (0, 1, -1)
Algorithmic connection for GKM graph with 3 nodes and valency 2
```

# Useful functions
```@docs
core
graph
num_edges
num_vertices
vertices_structure
GKMtools.lattice
GKMtools.flags
edges
vertices
label
degree
find_vertex_index
compact_flags
weight
other_vertex
print_labels
number_vertex_of_label
```