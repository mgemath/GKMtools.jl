@doc raw"""
    get_connection(gkm::AbstractGKM_graph) -> Union{Nothing, GKM_connection}

Return the connection of the given GKM graph if it is 3-independent, or if it is 2-valent and 2-independent,
 or if it has been set manually.
If one of the first two cases holds and the connection hasn't been calculated before,
it is saved in the `gkm` object for later use.
If none of the three cases hold, return `nothing`.

!!! note
    For GKM graphs of valency at least 3 that are not 3-independent, a connection may still exist,
    although uniqueness is not guaranteed.
    Use `get_any_connection` to get any compatible connection if one exists.

# Example
The unique connection for $\mathbb{P}^n$ has $\nabla_{(p\rightarrow q)}(p\rightarrow r)=(q\rightarrow r)$ for every triple of distinct vertices $(p, q, r)$, and $\nabla_{(p\rightarrow q)}(p\rightarrow q)=(q\rightarrow p)$ for every distinct vertices $p$ and $q$.

```julia-repl
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> C = get_connection(G)
GKM connection for GKM graph with 3 nodes and valency 2:
Connection:
Edge(1, 2) => [1, 2]
Edge(3, 1) => [2, 1]
Edge(2, 1) => [1, 2]
Edge(1, 3) => [2, 1]
Edge(3, 2) => [1, 2]
Edge(2, 3) => [1, 2]
a_i's:
Edge(1, 2) => ZZRingElem[2, 1]
Edge(3, 1) => ZZRingElem[2, 1]
Edge(2, 1) => ZZRingElem[2, 1]
Edge(1, 3) => ZZRingElem[1, 2]
Edge(3, 2) => ZZRingElem[1, 2]
Edge(2, 3) => ZZRingElem[1, 2]
```
"""
function get_connection(gkm::AbstractGKM_graph)::Union{Nothing, GKM_connection}
  if isnothing(gkm.connection)
    if (valency(gkm) >= 3 && is3_indep(gkm)) || (valency(gkm) == 2 && is2_indep(gkm)) || (valency(gkm)==1)
      gkm.connection = _build_GKM_connection(gkm)
    end
  end
  return gkm.connection
end

@doc raw"""
    get_any_connection(gkm::AbstractGKM_graph)::Union{Nothing, GKM_connection}

Return any connection for the given GKM graph, if there exists one, or `nothing` otherwise.
This connection is not guaranteed to have any special properties.
In particular, if `gkm` is the GKM graph of a sufficiently nice space, the returned connection
is not guaranteed to be the one induced by the geometry of the space.
"""
function get_any_connection(gkm::AbstractGKM_graph)::Union{Nothing, GKM_connection}
  con = get_connection(gkm)
  if !isnothing(con)
    return con
  elseif isnothing(gkm.anyConnection)
    gkm.anyConnection = _build_any_GKM_connection(gkm)
  end
  return gkm.anyConnection
end


@doc raw"""
    set_connection!(gkm::AbstractGKM_graph, con::GKM_connection)

Manually set the GKM connection of `gkm` to `con`.
This will overwrite any previously set connection.

# Example
After building the `GKM_connection` using `build_GKM_connection` like in the example above, we may assign it to the GKM graph using `set_connection!`:
```julia-repl
julia> G = projective_space(GKM_graph, 1);

julia> a = Dict{Edge, Vector{ZZRingElem}}();

julia> a[Edge(1, 2)] = [2];

julia> a[Edge(2, 1)] = [2];

julia> C = build_GKM_connection(G, a);

julia> set_connection!(G, C)
GKM connection for GKM graph with 2 nodes and valency 1:
Connection:
Edge(1, 2) => [1]
Edge(2, 1) => [1]
a_i's:
Edge(1, 2) => ZZRingElem[2]
Edge(2, 1) => ZZRingElem[2]
```
!!! note
    In this example, it is unnecessary to set the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function set_connection!(gkm::AbstractGKM_graph, con::GKM_connection)
  @req gkm == con.gkm "Connection belongs to the wrong GKM graph!"
  @req isvalid(con) "GKM connection is invalid!"
  
  gkm.connection = con
end

################
# Return the freshly calculated GKM_connection of the given GKM graph if it is unique.

# Warning:
#   1. This does not save the newly calculated GKM connection in the gkm object.
#   2. If the connection is unique or was set before, one should instead use get_connection().
#################
function _build_GKM_connection(gkm::AbstractGKM_graph) :: GKM_connection

  if valency(gkm) >= 3
    @req is3_indep(gkm) "GKM graph has valency >= 3 is not 3-independent"
  elseif valency(gkm) == 2
    @req is2_indep(gkm) "GKM graph has valency 2 and is not 2-independent"
  end

  @req isvalid(gkm) "GKM graph needs to be valid to build connection."

  # con[e][i] = j: along edge e from v to w, flag i at v connects to flag j at w
  con = Dict{Edge, Vector{Int64}}()
  val = valency(gkm)

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    con[e] = Vector{Int64}(undef, val)
    con[reverse(e)] = Vector{Int64}(undef, val)

    # get all flags at src(e)
    for i in 1:val
      wi = gkm.weights_at_vertex[s1][i]

      found = false
      # get all flags at dst(e)
      for j in 1:val
        wj = gkm.weights_at_vertex[s2][j]

        wdif = wi - wj # this will be a_i * w[e]

        if rank(matrix([ wdif; eW ])) == 1 # if true, flag j belongs to flag i.

          con[e][i] = j
          con[reverse(e)][j] = i
          found = true
          break
        end
      end
      @req found "No connection image found for edge $e, flag $(i)! The GKM graph does not admit a connection." # Assuming smoothness!
    end
  end

  res = build_GKM_connection(gkm, con)
  @req isvalid(res) "_build_GKM_connection built an invalid connection."
  return res
end

function _build_any_GKM_connection(gkm::AbstractGKM_graph) :: Union{Nothing, GKM_connection}

  @req isvalid(gkm) "GKM graph needs to be valid to build connection."

  if !is2_indep(gkm)
    println("Warning: The given GKM graph is not 2-independent!")
  end

  # con[e][i] = j: along edge e from v to w, flag i at v connects to flag j at w
  con = Dict{Edge, Vector{Int64}}()
  val = valency(gkm)

  # iterate over all unoriented edges
  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    con[e] = Vector{Int64}(undef, val)
    con[reverse(e)] = Vector{Int64}(undef, val)

    # make sure not to allocate some flag j to more than one flag i.
    allocatedFlags = Vector{Int64}()

    # get all flags at src(e)
    for i in 1:val
      wi = gkm.weights_at_vertex[s1][i]

      found = false
      # get all flags at dst(e)
      for j in 1:val

        (j in allocatedFlags) && continue

        wj = gkm.weights_at_vertex[s2][j]
        wdif = wi - wj # this will be a_i * w[e]

        if rank(matrix([ wdif; eW ])) == 1

          # Flag j is only a candidate for flag i if the resulting ai is an integer.
          aiIntegral::Bool = false
          for k in 1:rank(gkm.M)
            if eW[k] != 0
              tmp = wdif[k] // eW[k]
              aiIntegral = denominator(tmp) == 1
              break
            end
          end
          !aiIntegral && continue

          con[e][i] = j
          con[reverse(e)][j] = i
          push!(allocatedFlags, j)
          found = true
          break
        end
      end
      if !found
        println("No connection image found for edge $e, flag $(i)! The GKM graph does not admit a connection.")
        return nothing
      end
    end
  end
  c = build_GKM_connection(gkm, con)
  # If the below throws an error:
  # Is it because (e,e) !-> reverse(e)? This might happen only in the not 2-independent case.
  @req isvalid(c) "_build_any_GKM_connection built an invalid connection!"
  return c
end

@doc raw"""
    build_GKM_connection(gkm::AbstractGKM_graph, con::Dict{Edge, Vector{Int64}}) -> GKM_connection

Return the `GKM_connection` object (including information of the integers $a$) defined by the given connection map.
The connection map `con[e][i] = j` means that along edge `e` from `v`` to `w`, flag `i` at `v` connects to flag `j` at `w`.

!!! warning
    1. This function does not check whether the given connection map is valid (use `isvalid(::GKM_connection)` for that).
    2. This does not save the new connection to the gkm object (use `set_connection!(::AbstractGKM_graph, ::GKM_connection)` for that).

# Example
```julia-repl build_GKM_connection_from_a
julia> G = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> con = Dict{Edge, Vector{Int64}}()
Dict{Edge, Vector{Int64}}()

julia> con[Edge(1, 2)] = [1]
1-element Vector{Int64}: 1

julia> con[Edge(2, 1)] = [1]
1-element Vector{Int64}: 1

julia> C = build_GKM_connection(G, con)
GKM connection for GKM graph with 2 nodes and valency 1:
Connection:
Edge(1, 2) => [1]
Edge(2, 1) => [1]
a_i's:
Edge(1, 2) => ZZRingElem[2]
Edge(2, 1) => ZZRingElem[2]
```
!!! note
    In this example, it is unnecessary to define the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function build_GKM_connection(gkm::AbstractGKM_graph, con::Dict{Edge, Vector{Int64}}) :: GKM_connection
  a = connection_a_from_con(gkm, con)
  return GKM_connection(gkm, con, a)
end

@doc raw"""
    build_GKM_connection(gkm::AbstractGKM_graph, a::Dict{Edge, Vector{ZZRingElem}}) -> GKM_connection

Return the `GKM_connection` object (including the connection map $\nabla$) defined by the given integers `a`.
Here `a[e][i]` gives the integer coefficient such that `w[con[e][i] at dst(e)] = w[i at src(e)] - a[e][i] * w[e]`.

!!! warning
    1. This function does not check whether the given connection map is valid (use `isvalid(::GKM_connection)` for that).
    2. This does not save the new connection to the gkm object (use `set_connection!(::AbstractGKM_graph, ::GKM_connection)` for that).

# Example
```julia-repl build_GKM_connection_from_a
julia> G = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> a = Dict{Edge, Vector{ZZRingElem}}()
Dict{Edge, Vector{ZZRingElem}}()

julia> a[Edge(1, 2)] = [2]
1-element Vector{ZZRingElem}: 2

julia> a[Edge(2, 1)] = [2]
1-element Vector{ZZRingElem}: 2

julia> C = build_GKM_connection(G, a)
GKM connection for GKM graph with 2 nodes and valency 1:
Connection:
Edge(1, 2) => [1]
Edge(2, 1) => [1]
a_i's:
Edge(1, 2) => ZZRingElem[2]
Edge(2, 1) => ZZRingElem[2]
```
!!! note
    In this example, it is unnecessary to define the connection manually, since there is a unique one.
    To get it, simply use `get_connection(G)`.
"""
function build_GKM_connection(gkm::AbstractGKM_graph, a::Dict{Edge, Vector{ZZRingElem}}) :: GKM_connection
  con = connection_map_from_a(gkm, a)
  return GKM_connection(gkm, con, a)
end

# DEPRECATED: Compatibility function for old connection format
# Converts old Dict{Tuple{Edge, Edge}, ZZRingElem} to new Dict{Edge, Vector{ZZRingElem}} format
function build_GKM_connection(gkm::AbstractGKM_graph, a_old::Dict{Tuple{Edge, Edge}, ZZRingElem}) :: GKM_connection
  # Warning deactivated in order not to destroy doctests.
  # @warn "Using deprecated connection format Dict{Tuple{Edge, Edge}, ZZRingElem}. Please update to Dict{Edge, Vector{ZZRingElem}}." maxlog=1
  @req is_compact(gkm) "Deprecated build_GKM_connection is only valid for compact GKM spaces."

  # Convert old format to new format
  a_new = Dict{Edge, Vector{ZZRingElem}}()
  val = valency(gkm)

  for e in edges(gkm.g)
    a_new[e] = Vector{ZZRingElem}(undef, val)
    a_new[reverse(e)] = Vector{ZZRingElem}(undef, val)

    for i in 1:val
      # Find the corresponding entry in old format
      # In old format: a[(e, ei)] where ei is the i-th edge at src(e)
      edge_at_src = gkm.flag_to_edge[src(e)][i]
      @req !isnothing(edge_at_src) "Compact GKM graph wasn't compact after all... This should never happen!"
      a_new[e][i] = a_old[(e, edge_at_src)]

      edge_at_dst = gkm.flag_to_edge[dst(e)][i]
      @req !isnothing(edge_at_dst) "Compact GKM graph wasn't compact after all... This should never happen!"
      a_new[reverse(e)][i] = a_old[(reverse(e), edge_at_dst)]
    end
  end

  return build_GKM_connection(gkm, a_new)
end

# TODO: should probably remove the function below, as blowup connection needs to be redone properly.

# # Converter for blowup-specific connection format
# # In blowup code, connections are built as Dict{Tuple{Edge, Edge}, Edge}
# # where con[(e, flag_e)] = flag_e' means flag at src(e) represented by edge flag_e
# # connects to flag at dst(e) represented by edge flag_e'
# function _build_GKM_connection(gkm::AbstractGKM_graph, blowup_con::Dict{Tuple{Edge, Edge}, Edge}) :: GKM_connection
#   con_new = Dict{Edge, Vector{Int64}}()
#   val = valency(gkm)

#   for e in edges(gkm.g)
#     con_new[e] = Vector{Int64}(undef, val)
#     con_new[reverse(e)] = Vector{Int64}(undef, val)

#     for i in 1:val
#       # Find which edge represents flag i at src(e)
#       edge_at_i = gkm.flag_to_edge[src(e)][i]

#       if !isnothing(edge_at_i)
#         # Flag i is an edge flag
#         # Look up in blowup_con
#         if haskey(blowup_con, (e, edge_at_i))
#           # Get the destination edge
#           dst_edge = blowup_con[(e, edge_at_i)]
#           # Find which flag index at dst(e) this corresponds to
#           j = findfirst(edge -> edge == dst_edge, gkm.flag_to_edge[dst(e)])
#           if isnothing(j)
#             @error "Blowup connection inconsistency: destination edge not found in flags"
#             j = i  # fallback
#           end
#           con_new[e][i] = j
#         else
#           # No entry - use identity
#           con_new[e][i] = i
#         end
#       else
#         # Flag i is a standalone flag - use identity
#         con_new[e][i] = i
#       end

#       # Same for reverse direction
#       edge_at_i_rev = gkm.flag_to_edge[dst(e)][i]
#       if !isnothing(edge_at_i_rev)
#         if haskey(blowup_con, (reverse(e), edge_at_i_rev))
#           dst_edge = blowup_con[(reverse(e), edge_at_i_rev)]
#           j = findfirst(edge -> edge == dst_edge, gkm.flag_to_edge[src(e)])
#           if isnothing(j)
#             j = i
#           end
#           con_new[reverse(e)][i] = j
#         else
#           con_new[reverse(e)][i] = i
#         end
#       else
#         con_new[reverse(e)][i] = i
#       end
#     end
#   end

#   return build_GKM_connection(gkm, con_new)
# end

# Return the ai's belonging to the given GKM connection.
# Warning: This function does not check whether the given connection map is valid.
function connection_a_from_con(gkm::AbstractGKM_graph, con::Dict{Edge, Vector{Int64}}; check::Bool = true)::Dict{Edge, Vector{ZZRingElem}}

  a = Dict{Edge, Vector{ZZRingElem}}()
  val = valency(gkm)

  for e in edges(gkm.g)

    @req !is_zero(gkm.w[e]) "Weight zero edge found."

    s1 = src(e)
    s2 = dst(e)
    eW = gkm.w[e]

    a[e] = Vector{ZZRingElem}(undef, val)
    a[reverse(e)] = Vector{ZZRingElem}(undef, val)

    for i in 1:val

      j = con[e][i]
      wi = gkm.weights_at_vertex[s1][i]
      wj = gkm.weights_at_vertex[s2][j]
      wdif = wi - wj

      if check
        @req rank(matrix([ wdif; eW ])) == 1 "connection is incompatible with GKM graph"
      end

      # Find ai such that wdif = ai * eW
      # We need to find a non-zero component of eW
      ai::Union{Nothing, ZZRingElem} = nothing

      for k in 1:rank(gkm.M)
        if eW[k] != 0
          tmp = wdif[k] // eW[k]
          @req denominator(tmp) == 1 "GKM connection's a_i's must be integers!" # Assumption: x//y is integer if and only if denominator(x//y) == 1 in Oscar.
          ai = ZZ(tmp)
          break
        end
      end

      if isnothing(ai)
        # eW is zero, so wdif must also be zero (due to rank check)
        # In this case, ai is not well-defined - the connection a-value is arbitrary
        # This should not happen in a well-formed GKM graph
        error("Edge weight is zero for edge $e - cannot compute connection a-value")
      end

      a[e][i] = ai
      a[reverse(e)][j] = ai
    end
  end
  return a
end

# Build the connection map from the given collection of a's [cf. Liu--Sheshmani 2.(b) on p.4]
# Warning: The returned value is only unique if the GKM has no repeated weights at any vertex (which is required for it to be valid).
function connection_map_from_a(gkm::AbstractGKM_graph, a::Dict{Edge, Vector{ZZRingElem}})::Dict{Edge, Vector{Int64}}

  #TODO: check uniqueness as well!

  con = Dict{Edge, Vector{Int64}}()
  val = valency(gkm)

  for e in edges(gkm.g)

    s1 = src(e)
    s2 = dst(e)

    con[e] = Vector{Int64}(undef, val)
    con[reverse(e)] = Vector{Int64}(undef, val)

    for i in 1:val

      ai = a[e][i]
      wi = gkm.weights_at_vertex[s1][i]
      target_weight = wi - (gkm.weightType <: QQFieldElem ? QQ(ai) : ai) * gkm.w[e] # following [Liu--Sheshmani 2.(b) on p.4]

      resultFound = false

      for j in 1:val
        wj = gkm.weights_at_vertex[s2][j]

        if target_weight == wj

          con[e][i] = j
          con[reverse(e)][j] = i
          resultFound = true
          break
        end
      end
      @req resultFound "No flag found for edge $e, flag $i using connection a's"
    end
  end
  return con
end

@doc raw"""
    isvalid(con::GKM_connection; printDiagnostics::Bool=true) -> Bool

Return `true` if the given connection is valid for its GKM graph. This holds if and only if all of the following hold:
  1. `con.con[e]` and `con.a[e]` are set for all edges $e$ in the graph
  2. For each edge $e$, the flag corresponding to $e$ at src($e$) maps to the flag corresponding to reverse($e$) at dst($e$)
  3. For each edge $e$, the a-value for the flag corresponding to $e$ is 2
  4. Every flag connection satisfies the relation: `w[con[e][i] at dst(e)] = w[i at src(e)] - a[e][i] * w[e]`

# Example
```julia-repl isvalid_con
julia> G = projective_space(GKM_graph, 1);

julia> C = get_connection(G)
GKM connection for GKM graph with 2 nodes and valency 1

julia> C.con[Edge(1, 2)][1] = 2 # Invalid!
2

julia> isvalid(C)
Connection is invalid
false
```
"""
function isvalid(con::GKM_connection; printDiagnostics::Bool=true)::Bool

  if !isvalid(con.gkm)
    printDiagnostics && println("GKM graph of connection is invalid, so connection cannot be valid")
    return false
  end

  val = valency(con.gkm)

  for e_base in edges(con.gkm.g)
    for e in [e_base, reverse(e_base)]
      if !haskey(con.con, e)
        printDiagnostics && println("Connection misses key for e=$e.")
        return false
      elseif !haskey(con.a, e)
        printDiagnostics && println("Connection misses a for e=$e.")
        return false
      elseif length(con.con[e]) != val
        printDiagnostics && println("Connection for edge $e has wrong length: expected $val, got $(length(con.con[e])).")
        return false
      elseif length(con.a[e]) != val
        printDiagnostics && println("Connection a for edge $e has wrong length: expected $val, got $(length(con.a[e])).")
        return false
      end

      # Check that the flag corresponding to edge e maps correctly
      i_e = con.gkm.edge_to_flag_index[e]
      i_rev = con.gkm.edge_to_flag_index[reverse(e)]

      if con.con[e][i_e] != i_rev
        printDiagnostics && println("Connection doesn't map edge flag correctly: con[$e][$i_e] should be $i_rev, got $(con.con[e][i_e]).")
        return false
      end

      if con.a[e][i_e] != ZZ(2)
        printDiagnostics && println("Connection a-value for edge flag is not 2: a[$e][$i_e] = $(con.a[e][i_e]).")
        return false
      end
    end
  end

  # Check all flag connections satisfy the weight relation
  # Check both e and reverse(e) to ensure consistency
  for e_base in edges(con.gkm.g)
    for e in [e_base, reverse(e_base)]
      s1 = src(e)
      s2 = dst(e)
      eW = con.gkm.w[e]

      for i in 1:val
        j = con.con[e][i]
        ai = con.a[e][i]

        if j < 1 || j > val
          printDiagnostics && println("Connection maps flag $i at edge $e to invalid flag index $j.")
          return false
        end

        wi = con.gkm.weights_at_vertex[s1][i]
        wj = con.gkm.weights_at_vertex[s2][j]

        if wj != wi - base_ring(con.gkm.M)(ai) * eW
          printDiagnostics && println("Connection relation violated for edge $e, flag $i: w[$j at dst] != w[$i at src] - $ai * w[$e].")
          return false
        end
      end
    end
  end

  return true
end

function Base.show(io::IO, con::GKM_connection)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM connection")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", con::GKM_connection)

  print(io, "GKM connection for GKM graph with $(n_vertices(con.gkm.g)) nodes and valency $(valency(con.gkm)):")
  print(io, "\nConnection:")
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

function _get_connection_as(e::Edge, C::GKM_connection)
  G = C.gkm
  v = src(e)
  res = Vector{ZZRingElem}()
  for n in all_neighbors(G.g, v)
    push!(res, C.a[(e, Edge(v, n))])
  end
  return res
end

function _get_connection_as(src::String, dst::String, C::GKM_connection)
  e = edgeFromLabels(C.gkm, src, dst)
  return get_connection_as(e, C)
end