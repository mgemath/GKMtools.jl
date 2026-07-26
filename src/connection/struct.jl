struct Connection{R} <: AbstractGKMConnection{R}
  transport::Dict{Edge, Vector{Int}}
  a::Dict{Edge, Vector{R}}
  connection_type::String
end

connection(G::AbstractGKMGraph) = G.connection
empty_connection(R) = Connection{R}(Dict{Edge, Vector{Int}}(), Dict{Edge, Vector{R}}(), "")
is_empty(c::AbstractGKMConnection) = c.transport == Dict{Edge, Vector{Int}}()
# Interface

transport(G::AbstractGKMGraph) = transport(connection(G))
transport(c::AbstractGKMConnection) = c.transport

coefficients(G::AbstractGKMGraph) = coefficients(connection(G))
coefficients(c::AbstractGKMConnection) = c.a

connection_type(c::AbstractGKMConnection) = c.connection_type
connection_type(G::AbstractGKMGraph) = connection_type(connection(G))

is_canonical(c::AbstractGKMConnection) = !(is_empty(c.connection_type) || c.connection_type == "Algorithmic")
is_algorithmic(c::AbstractGKMConnection) = c.connection_type == "Algorithmic"


@doc"""
    print_connection(G::AbstractGKMGraph; verbose::Bool = true)

Print all relevant information of the connection of `G`.

# Examples
```jldoctest
julia> G = gkm_graph_of_toric(projective_space(NormalToricVariety, 3))
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 0, 0, 1)
3 -> 1 => (0, -1, 0, 1)
3 -> 2 => (1, -1, 0, 0)
4 -> 1 => (0, 0, -1, 1)
4 -> 2 => (1, 0, -1, 0)
4 -> 3 => (0, 1, -1, 0)
Algorithmic connection for GKM graph with 4 nodes and valency 3
```
Algorithmic connection means that the connection has been generated using the internal algorithm. 
This is perfectly fine since this GKM admits only one connection.
```jldoctest
julia> print_connection(G)
Algorithmic connection for GKM graph with 4 nodes and valency 3
Connection:
3 -> 4 => [1, 2, 3]
3 -> 1 => [1, 3, 2]
4 -> 1 => [1, 2, 3]
1 -> 3 => [1, 3, 2]
2 -> 1 => [2, 3, 1]
1 -> 2 => [3, 1, 2]
3 -> 2 => [1, 2, 3]
1 -> 4 => [1, 2, 3]
4 -> 2 => [2, 1, 3]
2 -> 3 => [1, 2, 3]
4 -> 3 => [1, 2, 3]
2 -> 4 => [2, 1, 3]
a_i's:
3 -> 4 => [1, 2, 1]
3 -> 1 => [1, 1, 2]
4 -> 1 => [1, 1, 2]
1 -> 3 => [1, 2, 1]
2 -> 1 => [1, 1, 2]
1 -> 2 => [2, 1, 1]
3 -> 2 => [2, 1, 1]
1 -> 4 => [1, 1, 2]
4 -> 2 => [2, 1, 1]
2 -> 3 => [2, 1, 1]
4 -> 3 => [1, 2, 1]
2 -> 4 => [1, 2, 1]
```
Here `3 -> 4 => [1, 2, 3]` means that transport along the edge from `3` to `4`
maps the `i`-th flag at vertex `3` to the `i`-th flag at vertex `4`. The entry
`3 -> 4 => [1, 2, 1]` under `a_i's` gives the corresponding coefficients:
for `i = 1, 2, 3`, the weight of the transported flag is the weight of the
`i`-th flag at vertex `3` minus, respectively, `1`, `2`, or `1` times the
weight of the edge from `3` to `4`.
"""
function print_connection(G::AbstractGKMGraph; verbose::Bool = true)
  print_connection(stdout, G; verbose)
end


function print_connection(io::IO, G::AbstractGKMGraph{R}; verbose::Bool = true) where R
  c = connection(G)
  if verbose
    print(io, "\n$(c.connection_type) connection for GKM graph with $(num_vertices(G)) nodes and valency $(valency(G))")
    print(io, "\nConnection:")
    for k in keys(c.transport)
      label_v = label(G, src(k))
      label_w = label(G, dst(k))
      print(io, "\n$label_v -> $label_w => $(transport(G)[k])")
      # print("\n$k => $(c.transport[k])")
    end
    print(io, "\na_i's:")
    for k in keys(c.a)
      label_v = label(G, src(k))
      label_w = label(G, dst(k))
      print(io, "\n$label_v -> $label_w => $(Int.(coefficients(G)[k]))")
      # print("\n$k => $(c.a[k])")
    end
    print(io, "\n")
  else
    if Oscar.is_terse(io)
      # no nested printing
      print(io, "\nGKM connection")
    else
      # nested printing allowed, preferably terse
      print(io, "\n$(c.connection_type) connection for GKM graph with $(num_vertices(G)) nodes and valency $(valency(G))")
    end
  end
end
