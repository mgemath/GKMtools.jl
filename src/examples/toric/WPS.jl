"""
    stacky_weighted_projective_space_fan(x::Vector{Int})

Takes a vector of positive integer weights `x` and returns:
1. `ray_mat`: An n x (n+1) matrix where the columns represent the rays.
2. `inc_mat`: A boolean Incidence matrix (rows = cones, columns = rays),
              ordered with bigger cones first.
"""
function stacky_weighted_projective_space_fan(x::Vector{Int})
    any(w -> w <= 0, x) && error("Weights must be positive integers.")
    isempty(x) && throw(ArgumentError("the weight vector must be non-empty"))

    n_plus_1 = length(x)
    n = n_plus_1 - 1

    # A simplicial fan for weighted projective space is obtained by taking the
    # standard basis rays in the first n coordinates and one extra ray given by
    # the negative of the remaining weights.
    ray_mat = zero_matrix(ZZ, n_plus_1, n)
    for i in 1:n
        ray_mat[i, i] = 1
    end

    last_ray = -Int.(x[2:end])
    divisor = gcd(abs.(last_ray)...)
    for j in 1:n
        ray_mat[n_plus_1, j] = div(last_ray[j], divisor)
    end

    # Incidence matrix: every maximal cone omits one ray, so it is indexed by
    # the subsets of rays of size k.
    num_cones = 2^n_plus_1 - 1
    inc_mat = zero!(BitArray(undef, (num_cones, n_plus_1)))
    row_idx = 1
    for k in n:-1:0
        combs = Combinatorics.combinations(1:n_plus_1, k)
        for comb in combs
            for col in comb
                inc_mat[row_idx, col] = true
            end
            row_idx += 1
        end
    end

    return WeightedProjectiveSpace(ray_mat, inc_mat, x)
end

function affine_space_stacky_cone(n::Int)
  n >= 0 || throw(ArgumentError("The dimension must be nonnegative."))

  return StackyCone(
    identity_matrix(ZZ, n),  # stacky rays e₁, …, eₙ
    identity_matrix(ZZ, n),  # primitive rays e₁, …, eₙ
    n,                       # rank of the ambient lattice
    collect(1:n)             # global ray indices
  )
end