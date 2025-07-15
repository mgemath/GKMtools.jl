@doc raw"""
    get_connection(G::GKM_graph; store::Bool = true) -> Union{Nothing, GKM_connection}

It returns the (unique) connection of the given GKM graph if it is 3-independent, or if it is 2-valent and 2-independent, or if it has been set manually.
If none of the three cases hold, return `nothing`. If `store` is true, the new connection is stored inside `G` for later use


!!! note
    For GKM graphs of valency at least 3 that are not 3-independent, a connection may still exist, although uniqueness is not guaranteed.
    Use `get_any_connection` to get any compatible connection if one exists.

# Example
The unique connection for $\mathbb{P}^n$ has $\nabla_{(p\rightarrow q)}(p\rightarrow r)=(q\rightarrow r)$ for every triple of distinct vertices $(p, q, r)$, and $\nabla_{(p\rightarrow q)}(p\rightarrow q)=(q\rightarrow p)$ for every distinct vertices $p$ and $q$.

```jldoctest
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> C = get_connection(G)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(3, 1), Edge(3, 2)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 3)) => Edge(2, 3)
(Edge(2, 3), Edge(2, 3)) => Edge(3, 2)
(Edge(2, 3), Edge(2, 1)) => Edge(3, 1)
(Edge(3, 2), Edge(3, 1)) => Edge(2, 1)
(Edge(1, 3), Edge(1, 2)) => Edge(3, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
(Edge(1, 3), Edge(1, 3)) => Edge(3, 1)
(Edge(2, 1), Edge(2, 3)) => Edge(1, 3)
(Edge(3, 2), Edge(3, 2)) => Edge(2, 3)
(Edge(3, 1), Edge(3, 1)) => Edge(1, 3)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(3, 1), Edge(3, 2)) => 1
(Edge(1, 2), Edge(1, 3)) => 1
(Edge(2, 3), Edge(2, 3)) => 2
(Edge(2, 3), Edge(2, 1)) => 1
(Edge(3, 2), Edge(3, 1)) => 1
(Edge(1, 3), Edge(1, 2)) => 1
(Edge(1, 2), Edge(1, 2)) => 2
(Edge(1, 3), Edge(1, 3)) => 2
(Edge(2, 1), Edge(2, 3)) => 1
(Edge(3, 2), Edge(3, 2)) => 2
(Edge(3, 1), Edge(3, 1)) => 2
```
"""
function get_connection(G::GKM_graph; store::Bool = true)::Union{Nothing, GKM_connection}
  
  if isnothing(G.connection)
    if (valency(G) >= 3 && is3_indep(G)) || (valency(G) == 2 && is2_indep(G)) || (valency(G)==1)
      connection = _build_gkm_connection(G)

      if store
        G.connection = connection
      end

    end
  end

  return G.connection
end

@doc raw"""
    get_any_connection(G::GKM_graph; store::Bool = true) -> Union{Nothing, GKM_connection}

It returns the connection of `G` if there exists one. Otherwise, it tries to construct one. This new connection is not guaranteed to have any special properties.
In particular, it is not guaranteed to be the one induced by the geometry of a space. If it is not possible to construct it, it returns `nothing`.
"""
function get_any_connection(G::GKM_graph; store::Bool = true)::Union{Nothing, GKM_connection}

  # try to the default connection
  con::Union{Nothing, GKM_connection} = get_connection(G; store)

  if !isnothing(con)
    return con # default connection exists
  end

  # default connection does not exists, try to construct one

  con = _build_any_gkm_connection(G)
  if store
    G.connection = con
  end

  return G.connection
end


@doc raw"""
    set_connection!(gkm::GKM_graph, con::GKM_connection)

Manually set the GKM connection of `gkm` to `con`.
This will overwrite any previously set connection.

# Example
After building the `GKM_connection` using `build_gkm_connection` like in the example above, we may assign it to the GKM graph using `set_connection!`:
```jldoctest
julia> G = projective_space(GKM_graph, 1);

julia> a = Dict{Tuple{Edge, Edge}, ZZRingElem}();

julia> a[(Edge(1, 2), Edge(1, 2))] = 2;

julia> a[(Edge(2, 1), Edge(2, 1))] = 2;

julia> C = build_gkm_connection(G, a);

julia> set_connection!(G, C)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(1, 2), Edge(1, 2)) => 2
```
!!! note
    In this example, it is unnecessary to set the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function set_connection!(G::GKM_graph, con::GKM_connection)
  # @req gkm == con.gkm "Connection belongs to the wrong GKM graph!"
  @req isvalid_connection(G, con) "GKM connection is invalid!"
  
  G.connection = con
end

################
# Return the freshly calculated GKM_connection of the given GKM graph if it is unique.

# Warning:
#   1. This does not save the newly calculated GKM connection in the gkm object.
#   2. If the connection is unique or was set before, one should instead use get_connection().
#################
function _build_gkm_connection(gkm::GKM_graph) :: GKM_connection
  
  if valency(gkm) >= 3
    @req is3_indep(gkm) "GKM graph has valency >= 3 is not 3-independent"
  elseif valency(gkm) == 2
    @req is2_indep(gkm) "GKM graph has valency 2 and is not 2-independent"
  end

  # assign to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  con = Dict{Tuple{Edge, Edge}, Edge}()

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    # get all edges at src(e)
    for v in all_neighbors(gkm.g, s1)
      ei = Edge(s1,v)

      # get all edges at dst(e), where epi stands for "e prime i"
      for w in all_neighbors(gkm.g, s2)
        epi = Edge(s2, w)

        wdif = gkm.w[ei] - gkm.w[epi] # this will be a_i * w[e]
        
        if isone(rank(matrix([ wdif; eW ]))) # if true, epi belongs to ei.

          con[(e, ei)] = epi
          con[(reverse(e), epi)] = ei
          break
        end
      end
      @req haskey(con, (e, ei)) "No connection image found for ($e, $ei)! The GKM graph does not admit a connection." # Assuming smoothness!
    end
  end

  return build_gkm_connection(gkm, con)
end

function _build_any_gkm_connection(gkm::GKM_graph) :: Union{Nothing, GKM_connection}

  if !is2_indep(gkm)
    println("Warning: The given GKM graph is not 2-independent!")
  end

  # assign to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  con = Dict{Tuple{Edge, Edge}, Edge}()

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    # make sure not to allocate some epi to more than one ei.
    allocatedEpis = Vector{Edge}()

    # get all edges at src(e)
    for v in all_neighbors(gkm.g, s1)
      ei = Edge(s1,v)

      # get all edges at dst(e), where epi stands for "e prime i"
      for w in all_neighbors(gkm.g, s2)
        epi = Edge(s2, w)

        wdif = gkm.w[ei] - gkm.w[epi] # this will be a_i * w[e]
        
        if isone(rank(matrix([ wdif; eW ]))) && !(epi in allocatedEpis)
          
          # Epi is only a candidate for ei if the resulting ai is an integer.
          aiIntegral::Bool = false
          for j in 1:rank(gkm.M)
            if eW[j] != 0
              tmp = wdif[j] // eW[j]
              aiIntegral = denominator(tmp) == 1
              break
            end
          end
          !aiIntegral && continue

          con[(e, ei)] = epi
          con[(reverse(e), epi)] = ei
          push!(allocatedEpis, epi)
          break
        end
      end
      if !haskey(con, (e, ei))
        println("No connection image found for ($e, $ei)! The GKM graph does not admit a connection.")
        return nothing
      end
    end
  end
  c = build_gkm_connection(gkm, con)
  # Of the below throws an error:
  # Is it because (e,e) !-> reverse(e)? This might happen only in the not 2-independent case.
  @req isvalid(c) "_build_any_gkm_connection build an invalid connection!"
  return c
end

@doc raw"""
    build_gkm_connection(gkm::GKM_graph, con::Dict{Tuple{Edge, Edge}, Edge}) -> GKM_connection

Return the `GKM_connection` object (including information of the integers $a$) defined by the given connection map.
!!! warning
    1. This function does not check whether the given connection map is valid (use `isvalid(::GKM_connection)` for that).
    2. This does not save the new connection to the gkm object (use `set_connection!(::GKM_graph, ::GKM_connection)` for that).

# Example
```jldoctest build_gkm_connection_from_a
julia> G = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> con = Dict{Tuple{Edge, Edge}, Edge}()
Dict{Tuple{Edge, Edge}, Edge}()

julia> con[(Edge(1, 2), Edge(1, 2))] = Edge(2, 1)
Edge(2, 1)

julia> con[(Edge(2, 1), Edge(2, 1))] = Edge(1, 2)
Edge(1, 2)

julia> C = build_gkm_connection(G, con)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(1, 2), Edge(1, 2)) => 2
```
!!! note
    In this example, it is unnecessary to define the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function build_gkm_connection(gkm::GKM_graph, con::Dict{Tuple{Edge, Edge}, Edge})::GKM_connection
  a = connection_a_from_con(gkm, con)
  return GKM_connection(con, a)
end

@doc raw"""
    build_gkm_connection(gkm::GKM_graph, a::Dict{Tuple{Edge, Edge}, ZZRingElem}) -> GKM_connection

Return the `GKM_connection` object (including the connection map $\nabla$) defined by the given integers `a`.

!!! warning
    1. This function does not check whether the given connection map is valid (use `isvalid(::GKM_connection)` for that).
    2. This does not save the new connection to the gkm object (use `set_connection!(::GKM_graph, ::GKM_connection)` for that).

# Example
```jldoctest build_gkm_connection_from_a
julia> G = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> a = Dict{Tuple{Edge, Edge}, ZZRingElem}()
Dict{Tuple{Edge, Edge}, ZZRingElem}()

julia> a[(Edge(1, 2), Edge(1, 2))] = 2
2

julia> a[(Edge(2, 1), Edge(2, 1))] = 2
2

julia> C = build_gkm_connection(G, a)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(1, 2), Edge(1, 2)) => 2
```
!!! note
    In this example, it is unnecessary to define the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function build_gkm_connection(gkm::GKM_graph, a::Dict{Tuple{Edge, Edge}, ZZRingElem}) :: GKM_connection
  con = connection_map_from_a(gkm, a)
  return GKM_connection(con, a)
end

# Return the ai's belonging to the given GKM connection.
# Warning: This function does not check whether the given connection map is valid.
function connection_a_from_con(gkm::GKM_graph, con::Dict{Tuple{Edge, Edge}, Edge}; check::Bool = true)::Dict{Tuple{Edge, Edge}, ZZRingElem}

  a = Dict{Tuple{Edge, Edge}, ZZRingElem}()
  for e in edges(gkm.g)
    
    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    for ei in [Edge(s1,v) for v in all_neighbors(gkm.g, s1)]
      
      epi = con[(e, ei)]
      wdif = gkm.w[ei] - gkm.w[epi]

      if check
        @req rank(matrix([ wdif; eW ])) == 1 "connection is incompatible with GKM graph"
      end

      ai::ZZRingElem = ZZ(0)

      for j in 1:rank(gkm.M)
        if eW[j] != 0
          tmp = wdif[j] // eW[j]
          @req denominator(tmp) == 1 "GKM connection's a_i's must be integers!" # Assumption: x//y is integer if and only if denominator(x//y) == 1 in Oscar.
          ai = ZZ(tmp)
          break
        end
      end

      a[(e, ei)] = ai
      a[(reverse(e), epi)] = ai
    end
  end
  return a
end

# Build the connection map from the given collection of a's [cf. Liu--Sheshmani 2.(b) on p.4]
# Warning: The returned value is only unique if the GKM has no repeated weights at any vertex (which is required for it to be valid).
function connection_map_from_a(gkm::GKM_graph, a::Dict{Tuple{Edge, Edge}, ZZRingElem})::Dict{Tuple{Edge, Edge}, Edge}

  #TODO: check uniqueness as well!

  con = Dict{Tuple{Edge, Edge}, Edge}()
  for e in edges(gkm.g)
    for ei in [Edge(src(e),v) for v in all_neighbors(gkm.g, src(e))]

      ai = a[(e, ei)]
      wEpi = gkm.w[ei] - ai * gkm.w[e] # following [Liu--Sheshmani 2.(b) on p.4]

      resultFound = false

      for epi in [Edge(dst(e),v) for v in all_neighbors(gkm.g, dst(e))]

        if wEpi == gkm.w[epi]

          con[(e, ei)] = epi
          con[(reverse(e), epi)] = ei
          resultFound = true
          break
        end
      end
      @req resultFound "No edge found for ($e,$ei) using connection a's"
    end
  end
  return con
end

@doc raw"""
    isvalid_connection(G::GKM_graph, con::GKM_connection; printDiagnostics::Bool=true)::Bool -> Bool

Return `true` if `con` is a valid connection of `G` is valid. This holds if and only if all of the followings hold:
  1. `con.con` and `con.a` are set for all `(Edge(v,w), Edge(v,u))` where `(v,w)` and `(v,u)` are edges in the graph.
  2. con maps every `(e,e)` to `reverse(e)`
  3. a maps every `(e,e)` to `2`
  4. Every pair of edges `(e,ei)` with same source satisfies the relation of the associated a's (see above), i.e. `con.gkm.w[ei'] = con.gkm.w[ei] - con.a[(e,ei)] * con.gkm.w[e]`

# Example
```jldoctest isvalid_con
julia> G = projective_space(GKM_graph, 1);

julia> C = get_connection(G)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(1, 2), Edge(1, 2)) => 2

julia> C.con[(Edge(1, 2), Edge(1, 2))] = Edge(1, 2) # Should be Edge(2, 1)!
Edge(1, 2)

julia> isvalid_connection(G, C)
Connection doesn't map (e,e) to reverse(e) for e=Edge(1, 2).
false
```
"""
function isvalid_connection(G::GKM_graph, con::GKM_connection; printDiagnostics::Bool=true)::Bool
  return _isvalid_connection(G, con; printDiagnostics)
end

@doc raw"""
    isvalid_connection(G::GKM_graph; printDiagnostics::Bool=true)::Bool -> Bool

Return `true` if the connection assigned to `G` is valid. This holds if and only if all of the followings hold:
  1. `con.con` and `con.a` are set for all `(Edge(v,w), Edge(v,u))` where `(v,w)` and `(v,u)` are edges in the graph.
  2. con maps every `(e,e)` to `reverse(e)`
  3. a maps every `(e,e)` to `2`
  4. Every pair of edges `(e,ei)` with same source satisfies the relation of the associated a's (see above), i.e. `con.gkm.w[ei'] = con.gkm.w[ei] - con.a[(e,ei)] * con.gkm.w[e]`

# Example
```jldoctest isvalid_con
julia> G = projective_space(GKM_graph, 1);

julia> C = get_connection(G)
Connection:
(Edge(2, 1), Edge(2, 1)) => Edge(1, 2)
(Edge(1, 2), Edge(1, 2)) => Edge(2, 1)
a_i's:
(Edge(2, 1), Edge(2, 1)) => 2
(Edge(1, 2), Edge(1, 2)) => 2

julia> C.con[(Edge(1, 2), Edge(1, 2))] = Edge(1, 2) # Should be Edge(2, 1)!
Edge(1, 2)

julia> isvalid_connection(G, C)
Connection doesn't map (e,e) to reverse(e) for e=Edge(1, 2).
false
```
"""
function isvalid_connection(G::GKM_graph; printDiagnostics::Bool=true)::Bool
  return _isvalid_connection(G, G.connection; printDiagnostics)
end

function _isvalid_connection(G::GKM_graph, con::Union{Nothing, GKM_connection}; printDiagnostics::Bool=true)::Bool

  con === nothing && return false

  for e in edges(G.g)
    if !haskey(con.con, (e,e))
      printDiagnostics && println("Connection misses key (e,e) for e=$e.")
      return false
    elseif !haskey(con.a, (e,e))
      printDiagnostics && println("Connection misses a(e,e) for e=$e.")
      return false
    elseif !haskey(con.con, (reverse(e),reverse(e)))
        printDiagnostics && println("Connection misses key (e,e) for e=$(reverse(e)).")
        return false
    elseif !haskey(con.a, (reverse(e),reverse(e)))
        printDiagnostics && println("Connection misses a(e,e) for e=$(reverse(e)).")
        return false
    elseif con.con[(e,e)] != reverse(e)
      printDiagnostics && println("Connection doesn't map (e,e) to reverse(e) for e=$e.")
      return false
    elseif con.con[(reverse(e),reverse(e))] != e
      printDiagnostics && println("Connection doesn't map (e,e) to reverse(e) for e=$(reverse(e)).")
      return false
    elseif con.a[(e,e)] != ZZ(2)
      printDiagnostics && println("Connection does not satisfy a(e,e)=2 for e=$e.")
      return false
    elseif con.a[(reverse(e), reverse(e))] != ZZ(2)
      printDiagnostics && println("Connection does not satisfy a(e,e)=2 for e=$(reverse(e)).")
      return false 
    end
  end

  for v in 1:n_vertices(G.g)
    for w in 1:n_vertices(G.g)
      (v == w) && continue
      e = Edge(v,w)
      if has_edge(G.g, e)
        for u in all_neighbors(G.g, v)

          ei = Edge(v, u)
          if !haskey(con.con, (e, ei))
            printDiagnostics && println("Connection map misses value for ($e, $ei).")
            return false
          elseif !haskey(con.a, (e, ei))
            printDiagnostics && println("Connection misses value for a($e, $ei).")
            return false
          end
          epi = con.con[(e,ei)]
          ai = con.a[(e,ei)]
          if G.w[epi] != G.w[ei] - base_ring(G.M)(ai) * G.w[e]
            printDiagnostics && println("Connection map and a(e,ei) is inconsistent for (e, ei)=($e, $ei).")
            return false
          end
        end
      end
    end
  end
  return true
end

function Base.show(io::IO, con::GKM_connection)

  print(io, "GKM connection")
  # if Oscar.is_terse(io)
  #   # no nested printing
  #   print(io, "GKM connection")
  # else
  #   # nested printing allowed, preferably terse
  #   print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm))")
  # end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", con::GKM_connection)

  # print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm)):")
  print(io, "Connection:")
  for k in keys(con.con)
    print(io, "\n$k => $(con.con[k])")
  end
  # show(io, MIME"text/plain"(), con.con)
  print(io, "\na_i's:")
  for k in keys(con.a)
    print(io, "\n$k => $(con.a[k])")
  end
  # show(io, MIME"text/plain"(), con.a)
  
end