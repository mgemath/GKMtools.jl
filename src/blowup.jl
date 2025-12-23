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
    The natural connection induced on the blowup is not yet implemented on this branch.
    However, this is not required for applications in Gromov-Witten theory.

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

Let us also blowup a non-compact GKM graph:
```jldoctest
julia> G = empty_gkm_graph(3, 2, ["1", "2", "3"]);

julia> g1, g2 = gens(G.M);

julia> add_edge!(G, 1, 2, g1);

julia> add_edge!(G, 2, 3, -g2);

julia> add_standalone_flag!(G, 1, g2);

julia> add_standalone_flag!(G, 1, -g1-g2);

julia> add_standalone_flag!(G, 2, g1+g2);

julia> add_standalone_flag!(G, 3, g1);

julia> add_standalone_flag!(G, 3, -g1-g2);

julia> set_connection!(G, get_any_connection(G)) # Set connection manually G is not 3-independent.

julia> S_cpct = gkm_subgraph_from_vertices(G, [1, 2])
GKM graph is valid but not 3-independent, so connections may not be unique.
GKM subgraph of:
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (0, 1)
Standalone flags:
1.2 => (0, 1)
1.3 => (-1, -1)
2.3 => (1, 1)
3.2 => (1, 0)
3.3 => (-1, -1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 0)

julia> B = blow_up(S)
GKM subgraph of:
GKM graph with 5 nodes, valency 3 and axial function:
[1>F3] -> [1>F2] => (1, 2)
[2>3] -> [1>F3] => (-1, 0)
[2>F3] -> [1>F2] => (-1, 0)
[2>F3] -> [2>3] => (-1, -2)
3 -> [2>3] => (0, 1)
Standalone flags:
[1>F2].2 => (0, 1)
[1>F3].2 => (-1, -1)
[2>F3].2 => (1, 1)
3.1 => (1, 0)
3.2 => (-1, -1)
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
[1>F3] -> [1>F2] => (1, 2)
[2>3] -> [1>F3] => (-1, 0)
[2>F3] -> [1>F2] => (-1, 0)
[2>F3] -> [2>3] => (-1, -2)

julia> # And another example where the subgraph that we blow up has standalone flags:

julia> S_non_cpct = gkm_subgraph_from_flags(G, [1, 2], [[1, 2], [1, 3]])
GKM subgraph of:
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (0, 1)
Standalone flags:
1.2 => (0, 1)
1.3 => (-1, -1)
2.3 => (1, 1)
3.2 => (1, 0)
3.3 => (-1, -1)
Subgraph:
GKM graph with 2 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
Standalone flags:
1.2 => (0, 1)
2.2 => (1, 1)

julia> B2 = blow_up(S_non_cpct)
GKM subgraph of:
GKM graph with 3 nodes, valency 3 and axial function:
[2>3] -> [1>F3] => (-1, 0)
3 -> [2>3] => (0, 1)
Standalone flags:
[1>F3].1 => (-1, -1)
[1>F3].2 => (0, 1)
[2>3].1 => (1, 1)
3.1 => (1, 0)
3.2 => (-1, -1)
Subgraph:
GKM graph with 2 nodes, valency 2 and axial function:
[2>3] -> [1>F3] => (-1, 0)
Standalone flags:
[1>F3].1 => (0, 1)
[2>3].1 => (1, 1)
```

!!! warning
    Blowups along non-closed subgraphs may result in invalid GKM graphs, like in the following example.
    The convention is that `Edge(s,d)` with `s` in the subgraph and `d` only creates an edge in the
    blowup if the flag for `Edge(s,d)` at `s` is not contained in the subgraph.

```jldoctest
julia> P3 = projective_space(GKM_graph, 3)
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> S = gkm_subgraph_from_flags(P3, [1, 2], [[1, 2], [1, 2]])
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Subgraph:
GKM graph with 2 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0, 0)
Standalone flags:
1.2 => (1, 0, -1, 0)
2.2 => (0, 1, -1, 0)

julia> B = blow_up(S)
The number of flags is not the same for all vertices
┌ Warning: Creating GKM subgraph of invalid gkm graph. This may result in undefined behavior.
└ @ GKMtools ~/julia_workspace/GKMtools.jl/src/GKMsubgraphs.jl:243
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
[2>4] -> [1>4] => (-1, 1, 0, 0)
4 -> [1>4] => (-1, 0, 0, 1)
4 -> [2>4] => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Standalone flags:
[1>4].1 => (1, 0, -1, 0)
[2>4].1 => (0, 1, -1, 0)
Subgraph:
GKM graph with 2 nodes, valency 2 and axial function:
[2>4] -> [1>4] => (-1, 1, 0, 0)
Standalone flags:
[1>4].1 => (1, 0, -1, 0)
[2>4].1 => (0, 1, -1, 0)
```
The warning above is thrown because the blowup GKM graph has non-constant valency, 
so creating the exceptional locus as GKM subgraph of an invalid GKM graph throws a warning.
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

  @req n != valency(sub) "Cannot blow up subgraph of equal valency as supergraph."

  # For each vertex in subgraph, find its normal flags (flags not in the subgraph)
  # A flag is "in the subgraph" if:
  # - It's an edge flag for an edge in the subgraph, OR
  # - It's a standalone flag that the subgraph explicitly contains

  normalFlags = Vector{Vector{Int64}}(undef, nvSub)  # For each v_sub: list of normal flag indices in super

  for v_sub in 1:nvSub
    normal = Int64[]

    for flag_idx in 1:n
      if !(flag_idx in gkmSub.flagDict[v_sub])
        push!(normal, flag_idx)
      end
    end

    normalFlags[v_sub] = normal
  end

  # Codimension, which is constant among subgraph vertices by preceding isvalid checks.
  c = length(normalFlags[1])

  # Build vertex mappings and labels
  labels = String[]
  exceptional_map = Dict{Tuple{Int64, Int64}, Int64}()  # (v_sub, normal_flag_idx) -> blowup vertex
  exceptional_flags = Dict{Int64, Vector{Int64}}() # exceptional_flags[v_blowup] -> flags at v_blowup in the exceptional locus
  non_sub_map = Dict{Int64, Int64}()  # v_super -> blowup vertex

  vertex_count = 0

  # Create exceptional vertices
  for v_sub in 1:nvSub
    v_super = vDict[v_sub]
    for flag_idx in normalFlags[v_sub]
      vertex_count += 1
      exceptional_map[(v_sub, flag_idx)] = vertex_count
      exceptional_flags[vertex_count] = Int64[]

      # Create label
      edge_at_flag = super.flag_to_edge[v_super][flag_idx]
      if !isnothing(edge_at_flag)
        @req src(edge_at_flag) == v_super "edge-flag relation broken"
        neighbor = dst(edge_at_flag)
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

        fi, fj = add_edge!(blowup, v_i, v_j, w_j - w_i)
        push!(exceptional_flags[v_i], fi)
        push!(exceptional_flags[v_j], fj)
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
    for i in 1:length(sub.weights_at_vertex[v_sub])
      if isnothing(sub.flag_to_edge[v_sub][i])
        # This is a standalone flag in the subgraph
        weight = sub.weights_at_vertex[v_sub][i]
        for j in 1:c
          v_exc = exceptional_map[(v_sub, nf[j])]
          fj = add_standalone_flag!(blowup, v_exc, weight)
          push!(exceptional_flags[v_exc], fj)
        end
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
        !isnothing(ef) && (src(ef) == s && dst(ef) == d)
      end, 1:n)
      @req !isnothing(flag_s) "flag for edge is missing"
      if flag_s in normalFlags[s_sub]
        add_edge!(blowup, exceptional_map[(s_sub, flag_s)], non_sub_map[d], w)
      end

    elseif !s_in_sub && d_in_sub
      # d in subgraph, s not
      d_sub = findfirst(==(d), vDict)
      flag_d = findfirst(i -> begin
        ef = super.flag_to_edge[d][i]
        !isnothing(ef) && (src(ef) == d && dst(ef) == s)
      end, 1:n)
      @req !isnothing(flag_d) "flag for edge is missing"
      if flag_d in normalFlags[d_sub]
        add_edge!(blowup, non_sub_map[s], exceptional_map[(d_sub, flag_d)], w)
      end

    else
      # Both in subgraph
      s_sub = findfirst(==(s), vDict)
      d_sub = findfirst(==(d), vDict)

      for flag_s in normalFlags[s_sub]
        flag_d = con.con[e][flag_s]
        @req flag_d in normalFlags[d_sub] "connection and subgraph were not compatible after all"
        vs = exceptional_map[(s_sub, flag_s)]
        vd = exceptional_map[(d_sub, flag_d)]
        fs, fd = add_edge!(blowup, vs, vd, w)
        push!(exceptional_flags[vs], fs)
        push!(exceptional_flags[vd], fd)
      end
    end
  end

  # Build connection: not yet implemented.

  # Exceptional divisor
  exceptional_vertices = [exceptional_map[(v_sub, f)] for v_sub in 1:nvSub for f in normalFlags[v_sub]]
  exceptional_flags = [exceptional_flags[v_exc] for v_exc in exceptional_vertices]
  gkmSubgraphBlowup = gkm_subgraph_from_flags(blowup, exceptional_vertices, exceptional_flags)

  return gkmSubgraphBlowup
end
