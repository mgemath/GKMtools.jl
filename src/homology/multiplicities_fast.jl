"""
    _multiplicities(H2::GKM_H2, edge_list::Vector{Edge}, beta::CurveClass)

Solve the positive integer edge-multiplicity equation directly when the
normalized `ray_sum` functional gives a finite total-degree bound. Fall back to
`_multiplicities_polymake` for degenerate homology data.

This file is intended to be included immediately after `curve_classes.jl`.
The old implementation in that file must first be renamed to
`_multiplicities_polymake`.
"""
function _multiplicities(H2::GKM_H2, edge_list::Vector{Edge}, beta::CurveClass)
  isempty(edge_list) && return Set{Vector{Int}}()

  h2_rank = rank(H2.H2)
  h2_rank == 0 && return _multiplicities_polymake(H2, edge_list, beta)
  length(H2.ray_sum) == h2_rank ||
    return _multiplicities_polymake(H2, edge_list, beta)

  edge_classes = Vector{Vector{ZZRingElem}}(undef, length(edge_list))
  for (i, edge) in enumerate(edge_list)
    edge_generator = gens(H2.edge_lattice)[H2.edge_to_gen[edge]]
    image = H2.quotient(edge_generator)
    edge_classes[i] = ZZ[image[j] for j in 1:h2_rank]
  end

  beta_coordinates = ZZ[beta[j] for j in 1:h2_rank]
  total_degree = sum(H2.ray_sum[j] * beta_coordinates[j] for j in 1:h2_rank)
  total_degree < 0 && return Set{Vector{Int}}()
  max_total = Int(floor(total_degree))
  length(edge_list) > max_total && return Set{Vector{Int}}()

  # The direct bound is valid only when ray_sum evaluates to at least one on
  # every column. This is guaranteed by the standard GKM_H2 constructor, but
  # the check preserves correctness for manually constructed homology data.
  for column in edge_classes
    evaluation = sum(H2.ray_sum[j] * column[j] for j in 1:h2_rank)
    evaluation >= 1 ||
      return _multiplicities_polymake(H2, edge_list, beta)
  end

  dual_rays = [QQ[ray[j] for j in 1:h2_rank] for ray in rays(H2.dual_cone)]
  result = Set{Vector{Int}}()
  current = zeros(Int, length(edge_list))
  dead_states = Set{Tuple{Int,Tuple{Vararg{ZZRingElem}},Int}}()

  _multiplicities_dfs!(
    result, current, edge_classes, beta_coordinates, dual_rays,
    1, max_total, dead_states,
  )
  return result
end

@inline function _multiplicity_dot(functional, class_coordinates)
  return sum(functional[j] * class_coordinates[j] for j in eachindex(class_coordinates))
end

@inline function _multiplicity_effective(remaining, dual_rays)
  return all(ray -> _multiplicity_dot(ray, remaining) >= 0, dual_rays)
end

function _multiplicities_dfs!(
  result::Set{Vector{Int}},
  current::Vector{Int},
  columns::Vector{Vector{ZZRingElem}},
  remaining::Vector{ZZRingElem},
  dual_rays,
  column_index::Int,
  remaining_degree::Int,
  dead_states::Set{Tuple{Int,Tuple{Vararg{ZZRingElem}},Int}},
)
  number_columns = length(columns)
  number_left = number_columns - column_index + 1
  remaining_degree < number_left && return false
  _multiplicity_effective(remaining, dual_rays) || return false

  state = (column_index, Tuple(remaining), remaining_degree)
  state in dead_states && return false

  column = columns[column_index]

  # The final variable is determined uniquely. Avoid an otherwise potentially
  # long integer loop and verify every homology coordinate exactly.
  if column_index == number_columns
    candidate = nothing
    for j in eachindex(remaining)
      coefficient = column[j]
      if iszero(coefficient)
        iszero(remaining[j]) || (push!(dead_states, state); return false)
        continue
      end
      q, r = divrem(remaining[j], coefficient)
      iszero(r) || (push!(dead_states, state); return false)
      if isnothing(candidate)
        candidate = q
      elseif candidate != q
        push!(dead_states, state)
        return false
      end
    end

    isnothing(candidate) && return _multiplicities_polymake_terminal!(
      result, current, columns, remaining, column_index, remaining_degree,
    )
    candidate_int = Int(candidate)
    if 1 <= candidate_int <= remaining_degree
      current[column_index] = candidate_int
      push!(result, copy(current))
      return true
    end
    push!(dead_states, state)
    return false
  end

  upper = remaining_degree - (number_left - 1)
  for ray in dual_rays
    denominator = _multiplicity_dot(ray, column)
    denominator > 0 || continue
    upper = min(upper, Int(floor(_multiplicity_dot(ray, remaining) / denominator)))
  end

  found = false
  next_remaining = similar(remaining)
  for multiplicity in 1:upper
    for j in eachindex(remaining)
      next_remaining[j] = remaining[j] - multiplicity * column[j]
    end
    _multiplicity_effective(next_remaining, dual_rays) || continue
    current[column_index] = multiplicity
    found |= _multiplicities_dfs!(
      result, current, columns, copy(next_remaining), dual_rays,
      column_index + 1, remaining_degree - multiplicity, dead_states,
    )
  end

  found || push!(dead_states, state)
  return found
end

# A zero final column is not expected for invariant curves in standard GKM_H2
# data. Returning no direct solution triggers no false positive; callers using
# unusual data should be handled by the validation/fallback above.
function _multiplicities_polymake_terminal!(result, current, columns, remaining,
                                             column_index, remaining_degree)
  return false
end
