struct Connection{R} <: AbstractGKMConnection{R}
  transport::Dict{Edge, Vector{Int}}
  a::Dict{Edge, Vector{R}}
  connection_type::String
end

connection(G::AbstractGKMGraph) = G.connection
empty_connection(R) = Connection{R}(Dict{Edge, Vector{Int}}(), Dict{Edge, Vector{R}}(), "")
is_empty(c::AbstractGKMConnection) = c.transport == Dict{Edge, Vector{Int}}()
# Interface

transport(c::AbstractGKMConnection) = c.transport
coefficients(c::AbstractGKMConnection) = c.a

connection_type(c::AbstractGKMConnection) = c.connection_type
connection_type(G::AbstractGKMGraph) = connection_type(connection(G))

is_canonical(c::AbstractGKMConnection) = !(is_empty(c.connection_type) || c.connection_type == "Algorithmic")
is_algorithmic(c::AbstractGKMConnection) = c.connection_type == "Algorithmic"

function print_connection(G::AbstractGKMGraph{R}; extended::Bool = true) where R
  c = connection(G)
  if extended
    print("\n$(c.connection_type) connection for GKM graph with $(num_vertices(G)) nodes and valency $(valency(G))")
    print("\nConnection:")
    for k in keys(c.transport)
      label_v = label(G, src(k))
      label_w = label(G, dst(k))
      print("\n$label_v -> $label_w => $(c.transport[k])")
      # print("\n$k => $(c.transport[k])")
    end
    print("\na_i's:")
    for k in keys(c.a)
      label_v = label(G, src(k))
      label_w = label(G, dst(k))
      print("\n$label_v -> $label_w => $(Int.(c.a[k]))")
      # print("\n$k => $(c.a[k])")
    end
  else
    if Oscar.is_terse(io)
      # no nested printing
      print("\nGKM connection")
    else
      # nested printing allowed, preferably terse
      print("\n$(c.connection_type) connection for GKM graph with $(num_vertices(G)) nodes and valency $(valency(G))")
    end
  end
end
