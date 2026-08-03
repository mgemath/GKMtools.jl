# Constructors

The main way to construct GKM graphs is using one of the constructors provided in the page [Examples](../Examples/Examples.md). The construction of a GKM graph from scratch is possible but not advised.

# Example
Here we construct the 2-dimensional projective space from scratch. 
```jldoctest scratch
julia> using GKMtools, Oscar
```
Coordinate torus T = (C*)³ acting on P²:

(λ₀, λ₁, λ₂) ⋅ [x₀:x₁:x₂]
    = [λ₀x₀ : λ₁x₁ : λ₂x₂].

The fixed points are
p₀ = [1:0:0], p₁ = [0:1:0], p₂ = [0:0:1].

Character lattice of T.
```jldoctest scratch
julia> M = free_module(ZZ, 3) # also possible to use QQ, but we want to see the integral structure
Free module of rank 3 over ZZ

julia> t₀, t₁, t₂ = gens(M)
3-element Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}:
 (1, 0, 0)
 (0, 1, 0)
 (0, 0, 1)

julia> characters = (t₀, t₁, t₂)
((1, 0, 0), (0, 1, 0), (0, 0, 1))
```
The one-dimensional T-orbits connect every pair of fixed points,
so the underlying graph is a triangle.
```jldoctest scratch
julia> g = Graph{Undirected}(3);

julia> add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 2, 3);
```
Vertex labels.

```jldoctest scratch
julia> vertex_labels = Vertex.([
           "p₀ = [1:0:0]",
           "p₁ = [0:1:0]",
           "p₂ = [0:0:1]",
       ])
3-element Vector{Vertex}:
 Vertex("p₀ = [1:0:0]")
 Vertex("p₁ = [0:1:0]")
 Vertex("p₂ = [0:0:1]")
```
At pᵢ, the tangent weight pointing toward pⱼ is tⱼ - tᵢ.
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
For each edge, record which flag at either endpoint belongs to it.
```jldoctest scratch
julia> edge_flags = Dict{Edge,Tuple{Int,Int}}();

julia> for e in edges(g)
           i = src(e)
           j = dst(e)

           α = characters[j] - characters[i]

           # Opposite orientations of the same invariant P¹ have
           # opposite tangent weights.
           push!(flags[i], FlagWeight{ZZRingElem}( α))
           push!(flags[j], FlagWeight{ZZRingElem}(-α))

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
  Vertex Labels: ["p₀ = [1:0:0]", "p₁ = [0:1:0]", "p₂ = [0:0:1]"]
  Flag weights defined at 3 vertices

```
Finally, construct the GKM graph.
```jldoctest scratch
julia> P2 = gkm_graph(core_data)
GKM graph with 3 nodes, valency 2 and axial function:
p₁ = [0:1:0] -> p₀ = [1:0:0] => (1, -1, 0)
p₂ = [0:0:1] -> p₀ = [1:0:0] => (1, 0, -1)
p₂ = [0:0:1] -> p₁ = [0:1:0] => (0, 1, -1)
Algorithmic connection for GKM graph with 3 nodes and valency 2
```

# Useful functions
```@docs
core
graph
num_edges
vertices_structure
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