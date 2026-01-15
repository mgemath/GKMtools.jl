###############################################################################
#
#   LaTeX/TikZ Drawing of GKM Graphs
#
###############################################################################

# Helper to format coordinates for TikZ output
_fmt(x::Float64) = string(round(x, digits=2))
_fmt(x::Int) = string(x)

@doc raw"""
    latex_drawing(G::AbstractGKM_graph; scale::Float64=2.0, vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false) -> String

Generate LaTeX/TikZ code that depicts the GKM graph `G`.

The vertex positions are computed automatically such that:
1. Vertices are placed at sufficient distance to avoid overlap.
2. Edges (or flags) whose weights are scalar multiples of each other are drawn parallel.

# Arguments
- `G::AbstractGKM_graph`: The GKM graph to draw.
- `scale::Float64=2.0`: Scaling factor for the drawing.
- `vertex_size::String="1.5pt"`: Size of vertex circles.
- `show_chern_numbers::Bool=false`: If true, label each edge with its Chern number and color edges green if Chern number ≤ valency(G) + 1.
- `show_vertex_labels::Bool=true`: If false, vertex labels are not displayed.
- `mark_irreducible::Bool=false`: If true, irreducible edges are drawn in red (takes precedence over green coloring).

# Returns
A string containing LaTeX/TikZ code. The code requires `\usepackage{tikz}` in the preamble.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 2);

julia> pos = [(0.0, 0.0), (2.0, 0.0), (1.0, 1.73)];

julia> println(latex_drawing(G, pos))
\\begin{tikzpicture}
% Vertices
\\node[draw, circle, fill=black, inner sep=1.5pt, label=left:{1}] (v1) at (0.0, 0.0) {};
\\node[draw, circle, fill=black, inner sep=1.5pt, label=right:{2}] (v2) at (2.0, 0.0) {};
\\node[draw, circle, fill=black, inner sep=1.5pt, label=above:{3}] (v3) at (1.0, 1.73) {};
% Edges
\\draw (v2) -- (v1);
\\draw (v3) -- (v1);
\\draw (v3) -- (v2);
\\end{tikzpicture}

julia> pos_int = [(0, 0), (2, 0), (1, 2)];

julia> println(latex_drawing(G, pos_int))
\\begin{tikzpicture}
% Vertices
\\node[draw, circle, fill=black, inner sep=1.5pt, label=left:{1}] (v1) at (0, 0) {};
\\node[draw, circle, fill=black, inner sep=1.5pt, label=right:{2}] (v2) at (2, 0) {};
\\node[draw, circle, fill=black, inner sep=1.5pt, label=above:{3}] (v3) at (1, 2) {};
% Edges
\\draw (v2) -- (v1);
\\draw (v3) -- (v1);
\\draw (v3) -- (v2);
\\end{tikzpicture}
```

The automatic layout can be obtained without positions:
```julia
julia> println(latex_drawing(G))  # automatic parallel-respecting layout
```
"""
function latex_drawing(G::AbstractGKM_graph; scale::Float64=2.0, vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false)::String
  positions = if has_attribute(G, :vert_pos)
    get_attribute(G, :vert_pos)
  else
    compute_gkm_layout(G, scale)
  end
  return generate_tikz_code(G, positions, vertex_size, show_chern_numbers, show_vertex_labels, mark_irreducible)
end

@doc raw"""
    latex_drawing(G::AbstractGKM_graph, positions::Vector{Tuple{Float64, Float64}}; vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false) -> String
    latex_drawing(G::AbstractGKM_graph, positions::Vector{Tuple{Int, Int}}; vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false) -> String

Generate LaTeX/TikZ code for the GKM graph `G` using user-provided vertex positions.

# Arguments
- `G::AbstractGKM_graph`: The GKM graph to draw.
- `positions`: A vector of (x, y) coordinates for each vertex. Can be `Vector{Tuple{Float64, Float64}}` or `Vector{Tuple{Int, Int}}`.
- `vertex_size::String="1.5pt"`: Size of vertex circles.
- `show_chern_numbers::Bool=false`: If true, label each edge with its Chern number and color edges green if Chern number ≤ valency(G) + 1.
- `show_vertex_labels::Bool=true`: If false, vertex labels are not displayed.
- `mark_irreducible::Bool=false`: If true, irreducible edges are drawn in red (takes precedence over green coloring).

# Example
```julia
G = projective_space(GKM_graph, 2)
pos = [(0, 0), (2, 0), (1, 2)]  # integer coordinates
println(latex_drawing(G, pos))

pos_float = [(0.0, 0.0), (2.0, 0.0), (1.0, 1.73)]  # float coordinates
println(latex_drawing(G, pos_float))

# With Chern numbers displayed
println(latex_drawing(G, pos, show_chern_numbers=true))
```
"""
function latex_drawing(G::AbstractGKM_graph, positions::Vector{Tuple{Float64, Float64}}; vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false)::String
  nv = n_vertices(G.g)
  @req length(positions) == nv "Number of positions must match number of vertices"
  return generate_tikz_code(G, positions, vertex_size, show_chern_numbers, show_vertex_labels, mark_irreducible)
end

function latex_drawing(G::AbstractGKM_graph, positions::Vector{Tuple{Int, Int}}; vertex_size::String="1.5pt", show_chern_numbers::Bool=false, show_vertex_labels::Bool=true, mark_irreducible::Bool=false)::String
  nv = n_vertices(G.g)
  @req length(positions) == nv "Number of positions must match number of vertices"
  return generate_tikz_code(G, positions, vertex_size, show_chern_numbers, show_vertex_labels, mark_irreducible)
end

"""
    generate_tikz_code(G::AbstractGKM_graph, positions::Vector{Tuple{T, T}}, vertex_size::String, show_chern_numbers::Bool, show_vertex_labels::Bool, mark_irreducible::Bool) -> String

Generate the TikZ code given vertex positions.

If `show_chern_numbers` is true, edges are labeled with their Chern numbers and edges
with Chern number at most `valency(G) + 1` are colored green.
If `mark_irreducible` is true, irreducible edges are colored red (takes precedence over green).
"""
function generate_tikz_code(G::AbstractGKM_graph, positions::Vector{Tuple{T, T}}, vertex_size::String, show_chern_numbers::Bool, show_vertex_labels::Bool, mark_irreducible::Bool)::String where T<:Real
  nv = n_vertices(G.g)
  val = valency(G)

  lines = String[]
  push!(lines, "\\begin{tikzpicture}")
  push!(lines, "% Vertices")

  # Determine label positions based on vertex positions
  for v in 1:nv
    x, y = positions[v]
    if show_vertex_labels
      label_pos = determine_label_position(positions, v)
      label = G.labels[v]
      push!(lines, "\\node[draw, circle, fill=black, inner sep=$vertex_size, label=$label_pos:{$label}] (v$v) at ($(_fmt(x)), $(_fmt(y))) {};")
    else
      push!(lines, "\\node[draw, circle, fill=black, inner sep=$vertex_size] (v$v) at ($(_fmt(x)), $(_fmt(y))) {};")
    end
  end

  push!(lines, "% Edges")

  # Compute label positions for edges to avoid overlapping labels
  edge_label_pos = if show_chern_numbers
    compute_edge_label_positions(G, positions)
  else
    Dict{Edge, Float64}()
  end

  for e in edges(G.g)
    s, d = src(e), dst(e)
    cc = curve_class(G, e)

    # Determine edge color: red for irreducible (highest priority), green for low Chern number, black otherwise
    is_irred = mark_irreducible && !is_reducible(G, cc)
    is_low_chern = show_chern_numbers && chern_number(G, cc) <= val + 1

    edge_style = if is_irred
      "red!70!black, thick"
    elseif is_low_chern
      "green!60!black, thick"
    else
      ""
    end

    if show_chern_numbers
      cn = chern_number(G, cc)
      pos_val = edge_label_pos[e]
      if isempty(edge_style)
        push!(lines, "\\draw (v$s) -- node[pos=$(_fmt(pos_val)), fill=white, inner sep=1pt] {\\small $cn} (v$d);")
      else
        push!(lines, "\\draw[$edge_style] (v$s) -- node[pos=$(_fmt(pos_val)), fill=white, inner sep=1pt] {\\small $cn} (v$d);")
      end
    else
      if isempty(edge_style)
        push!(lines, "\\draw (v$s) -- (v$d);")
      else
        push!(lines, "\\draw[$edge_style] (v$s) -- (v$d);")
      end
    end
  end

  # Draw standalone flags if any
  if !is_compact(G)
    push!(lines, "% Standalone flags")
    for v in 1:nv
      for (i, flag_edge) in enumerate(G.flag_to_edge[v])
        if isnothing(flag_edge)
          # Draw a short stub for standalone flag
          x, y = positions[v]
          weight = G.weights_at_vertex[v][i]
          # Compute direction from weight
          angle = weight_to_angle(weight)
          stub_length = 0.3
          x2 = x + stub_length * cos(angle)
          y2 = y + stub_length * sin(angle)
          push!(lines, "\\draw (v$v) -- ($(_fmt(x2)), $(_fmt(y2)));")
        end
      end
    end
  end

  push!(lines, "\\end{tikzpicture}")
  return join(lines, "\n")
end

"""
    determine_label_position(positions::Vector{Tuple{T, T}}, v::Int) -> String

Determine the best label position (above, below, left, right) for vertex v.
"""
function determine_label_position(positions::Vector{Tuple{T, T}}, v::Int)::String where T<:Real
  x, y = positions[v]

  # Compute centroid of all vertices
  cx = sum(p[1] for p in positions) / length(positions)
  cy = sum(p[2] for p in positions) / length(positions)

  # Place label away from center
  dx = x - cx
  dy = y - cy

  if abs(dy) > abs(dx)
    return dy > 0 ? "above" : "below"
  else
    return dx > 0 ? "right" : "left"
  end
end

"""
    compute_edge_label_positions(G::AbstractGKM_graph, positions::Vector{Tuple{T, T}}) -> Dict{Edge, Float64}

Compute label positions (as TikZ `pos` values between 0 and 1) for each edge.
Labels are shifted from the midpoint (0.5) if two edges would have overlapping labels.
"""
function compute_edge_label_positions(G::AbstractGKM_graph, positions::Vector{Tuple{T, T}})::Dict{Edge, Float64} where T<:Real
  edge_list = collect(edges(G.g))
  n_edges = length(edge_list)

  # Initialize all labels at midpoint
  label_pos = Dict{Edge, Float64}(e => 0.5 for e in edge_list)

  # Compute midpoints for all edges
  midpoints = Dict{Edge, Tuple{Float64, Float64}}()
  for e in edge_list
    s, d = src(e), dst(e)
    mx = (Float64(positions[s][1]) + Float64(positions[d][1])) / 2
    my = (Float64(positions[s][2]) + Float64(positions[d][2])) / 2
    midpoints[e] = (mx, my)
  end

  # Find edges with coinciding midpoints and shift their labels
  tolerance = 0.01
  for i in 1:n_edges
    for j in (i+1):n_edges
      e1, e2 = edge_list[i], edge_list[j]
      m1, m2 = midpoints[e1], midpoints[e2]

      # Check if midpoints coincide
      dist = sqrt((m1[1] - m2[1])^2 + (m1[2] - m2[2])^2)
      if dist < tolerance
        # Shift labels in opposite directions along their respective edges
        label_pos[e1] = 0.35
        label_pos[e2] = 0.65
      end
    end
  end

  return label_pos
end

"""
    weight_to_angle(weight) -> Float64

Convert a weight vector to an angle (in radians) for drawing direction.
Uses the first two components of the weight, or the first component and 0 if rank is 1.
"""
function weight_to_angle(weight)::Float64
  rk = rank(parent(weight))
  if rk == 0
    return 0.0
  elseif rk == 1
    w1 = weight[1]
    return w1 >= 0 ? 0.0 : π
  else
    w1, w2 = weight[1], weight[2]
    x = _to_float(w1)
    y = _to_float(w2)
    return atan(y, x)
  end
end

# Helper to convert ZZRingElem or QQFieldElem to Float64
function _to_float(val)::Float64
  if val isa QQFieldElem
    return Float64(numerator(val)) / Float64(max(1, denominator(val)))
  else
    return Float64(val)
  end
end

"""
    compute_weight_direction_groups(G::AbstractGKM_graph) -> Dict{Vector{Rational{BigInt}}, Vector{Edge}}

Group edges by their weight direction. Two weights are in the same group if one is a
scalar multiple of the other (i.e., they are parallel).

Returns a dictionary mapping normalized direction vectors to lists of edges.
"""
function compute_weight_direction_groups(G::AbstractGKM_graph)::Dict{Vector{Rational{BigInt}}, Vector{Edge}}
  groups = Dict{Vector{Rational{BigInt}}, Vector{Edge}}()

  for e in edges(G.g)
    weight = G.w[e]
    dir = normalize_weight_direction(weight)
    if !haskey(groups, dir)
      groups[dir] = Edge[]
    end
    push!(groups[dir], e)
  end

  return groups
end

"""
    normalize_weight_direction(weight) -> Vector{Rational{BigInt}}

Normalize a weight vector to a canonical direction representation.
Two weights that are scalar multiples of each other will have the same normalized direction.
"""
function normalize_weight_direction(weight)::Vector{Rational{BigInt}}
  rk = rank(parent(weight))
  if rk == 0
    return Rational{BigInt}[]
  end

  # Convert to rationals
  coeffs = Rational{BigInt}[]
  for i in 1:rk
    val = weight[i]
    if val isa QQFieldElem
      push!(coeffs, Rational{BigInt}(BigInt(numerator(val)), BigInt(denominator(val))))
    else
      push!(coeffs, Rational{BigInt}(BigInt(val)))
    end
  end

  # Find first non-zero coefficient
  first_nonzero = findfirst(c -> c != 0, coeffs)
  if isnothing(first_nonzero)
    return coeffs
  end

  # Normalize: divide by first non-zero, ensure it's positive
  factor = coeffs[first_nonzero]
  if factor < 0
    factor = -factor
    coeffs = [-c for c in coeffs]
  end

  return [c // abs(factor) for c in coeffs]
end

"""
    compute_gkm_layout(G::AbstractGKM_graph, scale::Float64) -> Vector{Tuple{Float64, Float64}}

Compute vertex positions for the GKM graph such that edges with parallel weights are drawn parallel.

This uses a constraint-based layout algorithm:
1. Group edges by weight direction (parallel weights = same group)
2. Assign each group a unique angle
3. Use force-directed placement to position vertices respecting the angle constraints
"""
function compute_gkm_layout(G::AbstractGKM_graph, scale::Float64)::Vector{Tuple{Float64, Float64}}
  nv = n_vertices(G.g)

  if nv == 0
    return Tuple{Float64, Float64}[]
  elseif nv == 1
    return [(0.0, 0.0)]
  elseif nv == 2
    return [(0.0, 0.0), (scale, 0.0)]
  end

  # Group edges by weight direction
  direction_groups = compute_weight_direction_groups(G)
  n_groups = length(direction_groups)

  if n_groups == 0
    # No edges: place vertices on a circle
    return circular_layout(nv, scale)
  end

  # Build edge-to-angle mapping
  # When rank_torus == 2, use the actual weight to determine the angle (slope proportional to weight)
  # Otherwise, spread angles evenly over [0, π)
  edge_angles = Dict{Edge, Float64}()
  rk = rank_torus(G)

  if rk == 2
    # Use weight vector directly to determine angle
    for (dir, edge_list) in direction_groups
      # dir is normalized, use it to compute angle
      angle = atan(Float64(dir[2]), Float64(dir[1]))
      for e in edge_list
        edge_angles[e] = angle
      end
    end
  else
    # For other ranks, spread angles evenly
    group_angles = Dict{Vector{Rational{BigInt}}, Float64}()
    for (i, dir) in enumerate(keys(direction_groups))
      group_angles[dir] = (i - 1) * π / n_groups
    end
    for (dir, edge_list) in direction_groups
      angle = group_angles[dir]
      for e in edge_list
        edge_angles[e] = angle
      end
    end
  end

  # Use constrained force-directed layout
  positions = constrained_force_layout(G, edge_angles, scale)

  return positions
end

"""
    circular_layout(n::Int, scale::Float64) -> Vector{Tuple{Float64, Float64}}

Place n vertices evenly on a circle.
"""
function circular_layout(n::Int, scale::Float64)::Vector{Tuple{Float64, Float64}}
  positions = Tuple{Float64, Float64}[]
  for i in 1:n
    angle = 2π * (i - 1) / n
    x = scale * cos(angle)
    y = scale * sin(angle)
    push!(positions, (x, y))
  end
  return positions
end

"""
    constrained_force_layout(G::AbstractGKM_graph, edge_angles::Dict{Edge, Float64}, scale::Float64) -> Vector{Tuple{Float64, Float64}}

Compute vertex positions using a force-directed algorithm with angle constraints.

Each edge must have a specific angle. The algorithm iteratively adjusts positions
to satisfy these constraints while also trying to:
- Keep vertices at reasonable distances
- Avoid vertex overlaps
"""
function constrained_force_layout(G::AbstractGKM_graph, edge_angles::Dict{Edge, Float64}, scale::Float64)::Vector{Tuple{Float64, Float64}}
  nv = n_vertices(G.g)

  # Initialize with circular layout
  positions = circular_layout(nv, scale)

  # Convert to mutable arrays
  x = [p[1] for p in positions]
  y = [p[2] for p in positions]

  # Iterative constraint solving
  n_iterations = 500
  learning_rate = 0.1
  min_dist = scale * 0.5  # Minimum distance between vertices

  for iter in 1:n_iterations
    # Compute forces
    fx = zeros(nv)
    fy = zeros(nv)

    # Repulsion between all vertex pairs
    for i in 1:nv
      for j in (i+1):nv
        dx = x[j] - x[i]
        dy = y[j] - y[i]
        dist = sqrt(dx^2 + dy^2)
        if dist < 0.01
          dist = 0.01
          dx, dy = rand() - 0.5, rand() - 0.5
        end

        # Repulsive force
        if dist < min_dist
          force = (min_dist - dist) / dist
          fx[i] -= force * dx
          fy[i] -= force * dy
          fx[j] += force * dx
          fy[j] += force * dy
        end
      end
    end

    # Apply forces
    for i in 1:nv
      x[i] += learning_rate * fx[i]
      y[i] += learning_rate * fy[i]
    end

    # Project edges onto their required angles
    for e in edges(G.g)
      s, d = src(e), dst(e)
      angle = edge_angles[e]

      # Current edge vector
      dx = x[d] - x[s]
      dy = y[d] - y[s]
      current_length = sqrt(dx^2 + dy^2)
      if current_length < 0.01
        current_length = scale * 0.5
      end

      # Target positions along the required angle
      # Keep midpoint fixed, adjust endpoints
      mx, my = (x[s] + x[d]) / 2, (y[s] + y[d]) / 2
      half_len = current_length / 2

      # Direction unit vector for required angle
      ux, uy = cos(angle), sin(angle)

      # Check which orientation is closer to current
      dot = dx * ux + dy * uy
      if dot < 0
        ux, uy = -ux, -uy
      end

      # New positions
      new_xs = mx - half_len * ux
      new_ys = my - half_len * uy
      new_xd = mx + half_len * ux
      new_yd = my + half_len * uy

      # Blend towards target
      blend = 0.3
      x[s] = (1 - blend) * x[s] + blend * new_xs
      y[s] = (1 - blend) * y[s] + blend * new_ys
      x[d] = (1 - blend) * x[d] + blend * new_xd
      y[d] = (1 - blend) * y[d] + blend * new_yd
    end
  end

  # Final pass: strictly enforce angle constraints
  # Do multiple passes to handle conflicts between edges sharing vertices
  for _ in 1:50
    for e in edges(G.g)
      s, d = src(e), dst(e)
      angle = edge_angles[e]

      dx = x[d] - x[s]
      dy = y[d] - y[s]
      current_length = sqrt(dx^2 + dy^2)
      if current_length < 0.01
        current_length = scale * 0.5
      end

      mx, my = (x[s] + x[d]) / 2, (y[s] + y[d]) / 2
      half_len = current_length / 2

      ux, uy = cos(angle), sin(angle)
      dot = dx * ux + dy * uy
      if dot < 0
        ux, uy = -ux, -uy
      end

      x[s] = mx - half_len * ux
      y[s] = my - half_len * uy
      x[d] = mx + half_len * ux
      y[d] = my + half_len * uy
    end
  end

  # Center the layout
  cx = sum(x) / nv
  cy = sum(y) / nv
  x .-= cx
  y .-= cy

  return [(x[i], y[i]) for i in 1:nv]
end
