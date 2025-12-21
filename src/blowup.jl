@doc"""
    blow_up(gkmSub::AbstractGKM_subgraph) -> AbstractGKM_subgraph

Return the tuple (GKM graph of blowup, GKM graph of exceptional divisor)
from (GKM graph, GKM subgraph, connection on supergraph), where both are encoded as AbstractGKM_subgraph.
!!! note
    The GKM graph needs to have the connection field set. The returned blowup graph and subgraph
    will also have the connection field set, but not the curveClasses field.
    (It will be calculated automatically on demand via `GKM_second_homology()`).
    Mathematically, this follows [GZ01; section 2.2.1](@cite).

!!! warning
    This will build an undirected graph. Behaviour with directed graphs as input is not tested.

# Examples
```jldoctest
julia> G = projective_space(GKM_graph, 3) # 3-dimensional projective space
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> S = gkm_subgraph_from_vertices(G, [1, 2]) # we take the subgraph of two vertices, it corresponds to a line
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1, 0, 0)

julia> blowupSub = blow_up(S) # blowup of P3 along the line S
GKM subgraph of:
GKM graph with 6 nodes, valency 3 and axial function:
[1>4] -> [1>3] => (0, 0, -1, 1)
[2>3] -> [1>3] => (-1, 1, 0, 0)
[2>4] -> [1>4] => (-1, 1, 0, 0)
[2>4] -> [2>3] => (0, 0, -1, 1)
3 -> [1>3] => (-1, 0, 1, 0)
3 -> [2>3] => (0, -1, 1, 0)
4 -> [1>4] => (-1, 0, 0, 1)
4 -> [2>4] => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
[1>4] -> [1>3] => (0, 0, -1, 1)
[2>3] -> [1>3] => (-1, 1, 0, 0)
[2>4] -> [1>4] => (-1, 1, 0, 0)
[2>4] -> [2>3] => (0, 0, -1, 1)

julia> Spoint = gkm_subgraph_from_vertices(G, [1]) # we take the subgraph of one vertex that is an invariant point
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 1 nodes, valency 0 and axial function:

julia> blowupPt = blow_up(Spoint) # blowup of P3 at a point
GKM subgraph of:
GKM graph with 6 nodes, valency 3 and axial function:
[1>3] -> [1>2] => (0, -1, 1, 0)
[1>4] -> [1>2] => (0, -1, 0, 1)
[1>4] -> [1>3] => (0, 0, -1, 1)
2 -> [1>2] => (-1, 1, 0, 0)
3 -> [1>3] => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> [1>4] => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 3 nodes, valency 2 and axial function:
[1>3] -> [1>2] => (0, -1, 1, 0)
[1>4] -> [1>2] => (0, -1, 0, 1)
[1>4] -> [1>3] => (0, 0, -1, 1)
```
"""
function Oscar.blow_up(gkmSub::AbstractGKM_subgraph)::AbstractGKM_subgraph

  # Get connection from supergraph (must use get_connection, not get_any_connection)
  con = get_connection(gkmSub.super)
  @req !isnothing(con) "Supergraph needs a connection from get_connection()"

  # Validate inputs
  @req isvalid(gkmSub) "invalid graph/subgraph pair"
  @req isvalid(con) "invalid connection"
  @req is_compatible_with_connection(gkmSub, con) "connection incompatible with subgraph"

  super = gkmSub.super
  sub = gkmSub.self
  nvSub = n_vertices(sub.g)
  nv = n_vertices(super.g)
  vDict = gkmSub.vDict
  M = super.M
  n = valency(super)  # constant valency

  # Trivial case: subgraph is entire graph
  if nvSub == nv
    return gkmSub
  end

  # For each vertex in subgraph, find its normal flags (flags not in the subgraph)
  # A flag is "in the subgraph" if:
  # - It's an edge flag for an edge in the subgraph, OR
  # - It's a standalone flag that the subgraph explicitly contains

  normalFlags = Vector{Vector{Int64}}(undef, nvSub)  # For each v_sub: list of normal flag indices in super

  for v_sub in 1:nvSub
    v_super = vDict[v_sub]
    normal = Int64[]

    for flag_idx in 1:n
      edge_at_flag = super.flag_to_edge[v_super][flag_idx]
      flag_in_subgraph = false

      if !isnothing(edge_at_flag)
        # Edge flag: in subgraph if the edge is in the subgraph
        neighbor = (src(edge_at_flag) == v_super) ? dst(edge_at_flag) : src(edge_at_flag)
        if has_edge(gkmSub, Edge(v_super, neighbor))
          flag_in_subgraph = true
        end
      else
        # Standalone flag: check if it's in the subgraph
        # A standalone flag at v_super is in the subgraph if sub contains it
        # We check by seeing if this flag's weight appears in sub.weights_at_vertex[v_sub]
        weight = super.weights_at_vertex[v_super][flag_idx]
        # Count how many standalone flags sub has at v_sub
        num_edges_at_v_sub = count(i -> !isnothing(sub.flag_to_edge[v_sub][i]), 1:length(sub.weights_at_vertex[v_sub]))
        # Check standalone flags in sub
        for i in (num_edges_at_v_sub + 1):length(sub.weights_at_vertex[v_sub])
          if sub.weights_at_vertex[v_sub][i] == weight
            flag_in_subgraph = true
            break
          end
        end
      end

      if !flag_in_subgraph
        push!(normal, flag_idx)
      end
    end

    normalFlags[v_sub] = normal
  end

  # Check constant codimension
  c = length(normalFlags[1])
  for v_sub in 2:nvSub
    @req length(normalFlags[v_sub]) == c "All vertices in subgraph must have the same number of normal flags (constant codimension)"
  end

  # Build vertex mappings and labels
  labels = String[]
  exceptional_map = Dict{Tuple{Int64, Int64}, Int64}()  # (v_sub, normal_flag_idx) -> blowup vertex
  non_sub_map = Dict{Int64, Int64}()  # v_super -> blowup vertex

  vertex_count = 0

  # Create exceptional vertices
  for v_sub in 1:nvSub
    v_super = vDict[v_sub]
    for flag_idx in normalFlags[v_sub]
      vertex_count += 1
      exceptional_map[(v_sub, flag_idx)] = vertex_count

      # Create label
      edge_at_flag = super.flag_to_edge[v_super][flag_idx]
      if !isnothing(edge_at_flag)
        neighbor = (src(edge_at_flag) == v_super) ? dst(edge_at_flag) : src(edge_at_flag)
        push!(labels, "[" * sub.labels[v_sub] * ">" * super.labels[neighbor] * "]")
      else
        push!(labels, "[" * sub.labels[v_sub] * ">F" * string(flag_idx) * "]")
      end
    end
  end

  # Non-subgraph vertices
  for v_super in 1:nv
    if !has_vertex(gkmSub, v_super)
      vertex_count += 1
      non_sub_map[v_super] = vertex_count
      push!(labels, super.labels[v_super])
    end
  end

  nvBlowup = vertex_count

  # Create blowup graph
  blowup = gkm_graph(Graph{Undirected}(nvBlowup), labels, M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{super.weightType}}(); checkLabels=false)

  # Add standalone flags to non-subgraph vertices
  for v_super in 1:nv
    if !has_vertex(gkmSub, v_super)
      v_blow = non_sub_map[v_super]
      # Add all standalone flags from this vertex in super
      for flag_idx in 1:n
        edge_at_flag = super.flag_to_edge[v_super][flag_idx]
        if isnothing(edge_at_flag)
          weight = super.weights_at_vertex[v_super][flag_idx]
          add_standalone_flag!(blowup, v_blow, weight)
        end
      end
    end
  end

  # Add edges

  # 1. Exceptional edges (complete graph at each blown-up vertex)
  for v_sub in 1:nvSub
    v_super = vDict[v_sub]
    nf = normalFlags[v_sub]

    for i in 1:c
      for j in (i+1):c
        v_i = exceptional_map[(v_sub, nf[i])]
        v_j = exceptional_map[(v_sub, nf[j])]

        w_i = super.weights_at_vertex[v_super][nf[i]]
        w_j = super.weights_at_vertex[v_super][nf[j]]

        add_edge!(blowup, v_i, v_j, w_j - w_i)
      end
    end

    # Add defining flag to each exceptional vertex as a standalone flag
    # (only if it was a standalone flag in the super graph)
    for i in 1:c
      v_exc = exceptional_map[(v_sub, nf[i])]
      flag_idx = nf[i]
      edge_at_flag = super.flag_to_edge[v_super][flag_idx]
      if isnothing(edge_at_flag)
        # This was a standalone flag in super, so add it as standalone in blowup
        defining_flag_weight = super.weights_at_vertex[v_super][flag_idx]
        add_standalone_flag!(blowup, v_exc, defining_flag_weight)
      end
      # If it was an edge flag, it will be added as an edge in the next section
    end

    # Add standalone flags from subgraph to all exceptional vertices
    # We need to check which standalone flags are actually in the subgraph
    num_edges_at_v_sub = count(i -> !isnothing(sub.flag_to_edge[v_sub][i]), 1:length(sub.weights_at_vertex[v_sub]))
    for i in (num_edges_at_v_sub + 1):length(sub.weights_at_vertex[v_sub])
      # This is a standalone flag in the subgraph
      weight = sub.weights_at_vertex[v_sub][i]
      for j in 1:c
        v_exc = exceptional_map[(v_sub, nf[j])]
        add_standalone_flag!(blowup, v_exc, weight)
      end
    end
  end

  # 2. Edges from original edges
  for e in edges(super.g)
    s = src(e)
    d = dst(e)
    w = super.w[e]

    s_in_sub = has_vertex(gkmSub, s)
    d_in_sub = has_vertex(gkmSub, d)

    if !s_in_sub && !d_in_sub
      # Both outside subgraph
      add_edge!(blowup, non_sub_map[s], non_sub_map[d], w)

    elseif s_in_sub && !d_in_sub
      # s in subgraph, d not
      s_sub = findfirst(==(s), vDict)
      # Find flag at s for this edge
      flag_s = findfirst(i -> begin
        ef = super.flag_to_edge[s][i]
        !isnothing(ef) && ((src(ef) == s && dst(ef) == d) || (dst(ef) == s && src(ef) == d))
      end, 1:n)
      if !isnothing(flag_s) && flag_s in normalFlags[s_sub]
        add_edge!(blowup, exceptional_map[(s_sub, flag_s)], non_sub_map[d], w)
      end

    elseif !s_in_sub && d_in_sub
      # d in subgraph, s not
      d_sub = findfirst(==(d), vDict)
      flag_d = findfirst(i -> begin
        ef = super.flag_to_edge[d][i]
        !isnothing(ef) && ((src(ef) == d && dst(ef) == s) || (dst(ef) == d && src(ef) == s))
      end, 1:n)
      if !isnothing(flag_d) && flag_d in normalFlags[d_sub]
        add_edge!(blowup, non_sub_map[s], exceptional_map[(d_sub, flag_d)], w)
      end

    else
      # Both in subgraph
      s_sub = findfirst(==(s), vDict)
      d_sub = findfirst(==(d), vDict)

      for flag_s in normalFlags[s_sub]
        flag_d = con.con[e][flag_s]
        if flag_d in normalFlags[d_sub]
          add_edge!(blowup, exceptional_map[(s_sub, flag_s)], exceptional_map[(d_sub, flag_d)], w)
        end
      end
    end
  end

  # Build connection
  # The connection must respect the actual flag ordering at each vertex
  # We need to build newCon based on ACTUAL edges, not assumed ordering

  # Try to use get_any_connection to compute a valid connection
  blowup_con = nothing
  try
    blowup_con = get_any_connection(blowup)
    if !isnothing(blowup_con)
      set_connection!(blowup, blowup_con)
    end
  catch e
    println("Warning: Could not build connection for blowup: $e")
  end

  # Exceptional divisor
  exceptional_vertices = [exceptional_map[(v_sub, f)] for v_sub in 1:nvSub for f in normalFlags[v_sub]]
  gkmSubgraphBlowup = gkm_subgraph_from_vertices(blowup, exceptional_vertices; include_standalone_flags=true)

  return gkmSubgraphBlowup
end
