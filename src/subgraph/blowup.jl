@doc raw"""
    blow_up(S::GKMSubgraph)
    blow_up(S::GKMSubgraph, weights::AbstractVector{<:Integer})

Blow up a smooth GKM graph along `S` and return the exceptional divisor as a
`GKMSubgraph` of the blowup. The weighted form assigns the positive integers
`weights` to the normal flags at the first vertex of `S`; the ambient
connection transports this assignment over the center. Weighted blowups preserve
the ambient character lattice.

# Examples
```jldoctest
julia> G = projective_space(GKMGraph, 3) # 3-dimensional projective space
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3

julia> S = subgraph_from_vertices(G, [1, 2]) # we take the subgraph of two vertices, it corresponds to a line
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1, 0, 0)
Restricted connection for GKM graph with 2 nodes and valency 1

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

julia> Spoint = subgraph_from_vertices(G, [1]) # we take the subgraph of one vertex that is an invariant point
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
function Oscar.blow_up(S::GKMSubgraph)
  normal = _normal_flags(S)
  isempty(normal) && throw(ArgumentError("the center has no normal directions"))
  return _blow_up(S, ones(Int, length(first(normal))))
end

function Oscar.blow_up(
  S::GKMSubgraph,
  weights::AbstractVector{<:Integer},
)
  return _blow_up(S, Int.(weights))
end

@doc raw"""
    weighted_blow_up(S::GKMSubgraph, weights::AbstractVector{<:Integer}) -> OrbifoldGKMSubgraph

See [`blow_up`](@ref).
"""
weighted_blow_up(S::GKMSubgraph, weights::AbstractVector{<:Integer}) =
  blow_up(S, weights)

function _normal_flags(S::GKMSubgraph)
  ambient = ambient_graph(S)
  return [
    setdiff(
      collect(eachindex(flags(ambient, vertex_to_ambient(S, v)))),
      S.flag_to_ambient[v],
    )
    for v in 1:num_vertices(subgraph(S))
  ]
end

function _blow_up(S::GKMSubgraph{R}, weights::Vector{Int}) where {R}
  all(>(0), weights) || throw(ArgumentError("blowup weights must be positive"))
  normal = _normal_flags(S)
  isempty(normal) && throw(ArgumentError("the center is empty"))
  codimension = length(first(normal))
  codimension > 0 || throw(ArgumentError("the center has no normal directions"))
  length(weights) == codimension || throw(DimensionMismatch(
    "expected $codimension blowup weights, got $(length(weights))",
  ))
  all(length(n) == codimension for n in normal) || throw(ArgumentError(
    "the center does not have constant codimension",
  ))

  transported_weights = _transport_normal_weights(S, normal, weights)
  # A non-uniform weighted exceptional fibre has characters
  # alpha_j / w_j - alpha_i / w_i.  They need not lie in the integral
  # character lattice even when the ambient graph does, so construct the
  # blowup over QQ in that case.
  C = all(isone, weights) ? R : QQFieldElem
  return _construct_blowup(S, normal, transported_weights, C)
end

function _transport_normal_weights(S, normal, initial_weights)
  center = subgraph(S)
  ambient = ambient_graph(S)
  con = connection(ambient)
  assigned = [Dict{Int,Int}() for _ in 1:num_vertices(center)]
  for (flag, multiplicity) in zip(normal[1], initial_weights)
    assigned[1][flag] = multiplicity
  end

  queue = [1]
  visited = falses(num_vertices(center))
  visited[1] = true
  while !isempty(queue)
    v = popfirst!(queue)
    for edge in edges(center)
      if src(edge) == v
        w, oriented = dst(edge), edge
      elseif dst(edge) == v
        w, oriented = src(edge), reverse(edge)
      else
        continue
      end
      ambient_edge = edge_to_ambient(S, oriented)
      for flag in normal[v]
        target = transport(con)[ambient_edge][flag]
        target in normal[w] || throw(ArgumentError(
          "the ambient connection does not preserve the normal bundle along $ambient_edge",
        ))
        multiplicity = assigned[v][flag]
        old = get(assigned[w], target, multiplicity)
        old == multiplicity || throw(ArgumentError(
          "blowup weights are not invariant under normal transport",
        ))
        assigned[w][target] = multiplicity
      end
      if !visited[w]
        visited[w] = true
        push!(queue, w)
      end
    end
  end
  all(visited) || throw(ArgumentError("the blowup center must be connected"))
  all(length(assigned[v]) == length(normal[v]) for v in eachindex(normal)) ||
    throw(ArgumentError("could not transport every normal weight"))
  return assigned
end

function _construct_blowup(S::GKMSubgraph, normal, multiplicities, ::Type{C}) where {C}
  ambient = ambient_graph(S)
  center = subgraph(S)
  old_lattice = lattice(ambient)
  M = free_module(parent(zero(C)), rank(old_lattice))
  basis = gens(M)
  lift_weight(w) = sum(i -> C(w[i]) * basis[i], eachindex(basis); init=zero(M))

  exceptional_map = Dict{Tuple{Int,Int},Int}()
  outside_map = zeros(Int, num_vertices(ambient))
  blowup_labels = BlowupVertex[]
  for v in 1:num_vertices(center)
    ambient_vertex = vertex_to_ambient(S, v)
    for flag in normal[v]
      exceptional_map[(v, flag)] = length(blowup_labels) + 1
      push!(blowup_labels, BlowupVertex(
        "[" * label(center, v) * ">" * string(flag) * "]",
        labels(ambient)[ambient_vertex],
      ))
    end
  end
  exceptional_count = length(blowup_labels)
  for v in vertices(ambient)
    has_ambient_vertex(S, v) && continue
    outside_map[v] = length(blowup_labels) + 1
    push!(blowup_labels, BlowupVertex(label(ambient, v), labels(ambient)[v]))
  end

  g = Graph{Undirected}(length(blowup_labels))
  new_flags = [FlagWeight{C}[] for _ in blowup_labels]
  edge_flags = Dict{Edge,Tuple{Int,Int}}()
  exceptional_flag_images = [Int[] for _ in 1:exceptional_count]

  function add_blowup_edge!(u, v, source_weight; exceptional_at_source=false, exceptional_at_target=false)
    add_edge!(g, u, v)
    push!(new_flags[u], FlagWeight{C}(source_weight))
    push!(new_flags[v], FlagWeight{C}(-source_weight))
    e = Edge(u, v)
    edge_flags[e] = (length(new_flags[u]), length(new_flags[v]))
    exceptional_at_source && push!(exceptional_flag_images[u], length(new_flags[u]))
    exceptional_at_target && push!(exceptional_flag_images[v], length(new_flags[v]))
    return e
  end

  # Exceptional weighted projective fibres.
  for v in 1:num_vertices(center)
    ambient_vertex = vertex_to_ambient(S, v)
    nf = normal[v]
    for i in 1:(length(nf) - 1), j in (i + 1):length(nf)
      fi, fj = nf[i], nf[j]
      if C === QQFieldElem
        wi = inv(C(multiplicities[v][fj])) * lift_weight(weight(ambient, ambient_vertex, fj))
        wj = inv(C(multiplicities[v][fi])) * lift_weight(weight(ambient, ambient_vertex, fi))
      else
        # Uniform weights differ from the ordinary projective weights by one
        # global scalar, so retain the integral model.
        wi = C(multiplicities[v][fi]) * lift_weight(weight(ambient, ambient_vertex, fj))
        wj = C(multiplicities[v][fj]) * lift_weight(weight(ambient, ambient_vertex, fi))
      end
      add_blowup_edge!(
        exceptional_map[(v, fi)], exceptional_map[(v, fj)], wi - wj;
        exceptional_at_source=true, exceptional_at_target=true,
      )
    end
  end

  # Proper transforms of ambient edges.
  for edge in edges(ambient)
    s, d = src(edge), dst(edge)
    local_s, local_d = ambient_to_vertex(S, s), ambient_to_vertex(S, d)
    s_in, d_in = !iszero(local_s), !iszero(local_d)
    source_flag, target_flag = _edge_flag_indices(ambient, edge)
    original_weight = lift_weight(weight(ambient, edge))

    if !s_in && !d_in
      add_blowup_edge!(outside_map[s], outside_map[d], original_weight)
    elseif s_in && !d_in
      source_flag in normal[local_s] || continue
      m = multiplicities[local_s][source_flag]
      add_blowup_edge!(
        exceptional_map[(local_s, source_flag)], outside_map[d], original_weight,
      )
    elseif !s_in && d_in
      target_flag in normal[local_d] || continue
      m = multiplicities[local_d][target_flag]
      add_blowup_edge!(
        outside_map[s], exceptional_map[(local_d, target_flag)], original_weight,
      )
    else
      ambient_edge = edge_to_ambient(S, Edge(local_s, local_d))
      for source_normal in normal[local_s]
        target_normal = transport(connection(ambient))[ambient_edge][source_normal]
        target_normal in normal[local_d] || throw(ArgumentError(
          "the center is not compatible with the ambient connection",
        ))
        add_blowup_edge!(
          exceptional_map[(local_s, source_normal)],
          exceptional_map[(local_d, target_normal)],
          original_weight;
          exceptional_at_source=true, exceptional_at_target=true,
        )
      end
    end
  end

  # Standalone flags are copied to unchanged vertices and lifted over the center.
  attached = [falses(length(flags(ambient, v))) for v in vertices(ambient)]
  for edge in edges(ambient)
    i, j = _edge_flag_indices(ambient, edge)
    attached[src(edge)][i] = true
    attached[dst(edge)][j] = true
  end
  for v in vertices(ambient)
    if iszero(ambient_to_vertex(S, v))
      for i in eachindex(attached[v])
        attached[v][i] && continue
        push!(new_flags[outside_map[v]], FlagWeight{C}(lift_weight(weight(ambient, v, i))))
      end
    else
      local_v = ambient_to_vertex(S, v)
      center_flags = Set(S.flag_to_ambient[local_v])
      for normal_flag in normal[local_v]
        exceptional_vertex = exceptional_map[(local_v, normal_flag)]
        for ambient_flag in center_flags
          attached[v][ambient_flag] && continue
          push!(new_flags[exceptional_vertex], FlagWeight{C}(lift_weight(weight(ambient, v, ambient_flag))))
          push!(exceptional_flag_images[exceptional_vertex], length(new_flags[exceptional_vertex]))
        end
        if !attached[v][normal_flag]
          m = multiplicities[local_v][normal_flag]
          normal_weight = lift_weight(weight(ambient, v, normal_flag))
          C === QQFieldElem && (normal_weight = inv(C(m)) * normal_weight)
          push!(new_flags[exceptional_vertex], FlagWeight{C}(
            normal_weight,
          ))
        end
      end
    end
  end

  typed_labels = Vector{BlowupVertex{eltype(GKMtools.labels(ambient))}}(blowup_labels)
  core_data = GKMCombinatorialData{C,eltype(typed_labels),FlagWeight{C}}(
    g, M, typed_labels, new_flags, edge_flags,
  )
  blowup_connection = build_gkm_connection(
    core_data;
    connection_type="Blowup connection of $(connection_type(ambient))",
  )
  blowup_graph = GKMGraph{C,eltype(typed_labels),FlagWeight{C}}(
    core_data,
    blowup_connection,
    create_cohomology(rank(M), length(typed_labels)),
    nothing,
  )

  exceptional_vertices = collect(1:exceptional_count)
  return _subgraph_from_selected_flags(
    blowup_graph, exceptional_vertices, exceptional_flag_images,
  )
end

function _subgraph_from_selected_flags(ambient::AbstractGKMGraph{C,V,F}, vertex_map, selected_flags) where {C,V,F}
  local_of = Dict(v => i for (i, v) in enumerate(vertex_map))
  g = Graph{Undirected}(length(vertex_map))
  local_flags = [F[flags(ambient, v)[i] for i in selected_flags[k]] for (k, v) in enumerate(vertex_map)]
  edge_flags = Dict{Edge,Tuple{Int,Int}}()

  for edge in edges(ambient)
    haskey(local_of, src(edge)) && haskey(local_of, dst(edge)) || continue
    ambient_source, ambient_target = _edge_flag_indices(ambient, edge)
    local_source = findfirst(==(ambient_source), selected_flags[local_of[src(edge)]])
    local_target = findfirst(==(ambient_target), selected_flags[local_of[dst(edge)]])
    (isnothing(local_source) || isnothing(local_target)) && continue
    s, d = local_of[src(edge)], local_of[dst(edge)]
    add_edge!(g, s, d)
    edge_flags[Edge(s, d)] = (local_source, local_target)
  end

  local_labels = labels(ambient)[vertex_map]
  data = GKMCombinatorialData{C,V,F}(
    g, lattice(ambient), local_labels, local_flags, edge_flags,
  )
  con = isempty(edge_flags) ? empty_connection(C) :
        build_gkm_connection(data; connection_type="Exceptional divisor")
  local_graph = GKMGraph{C,V,F}(
    data, con, create_cohomology(rank(lattice(ambient)), length(vertex_map)), nothing,
  )
  return GKMSubgraph(ambient, local_graph, vertex_map, selected_flags)
end
