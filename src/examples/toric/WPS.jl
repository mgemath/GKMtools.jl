"""
    stacky_weighted_projective_space_fan(x::Vector{Int})

Takes a vector of positive integer weights `x` and returns:
1. `ray_mat`: An n x (n+1) matrix where the columns represent the rays.
2. `inc_mat`: A boolean Incidence matrix (rows = cones, columns = rays),
              ordered with bigger cones first.
"""
function stacky_weighted_projective_space_fan(x::Vector{Int})
    any(w -> w <= 0, x) && error("Weights must be positive integers.")

    n_plus_1 = length(x)
    n = n_plus_1 - 1
    
    # 1. Natively allocate the exact output matrix: (n+1) rows by n columns.
    # Because Julia is column-major, updating the columns of this matrix is very fast.
    ray_mat = zero_matrix(ZZ, n_plus_1, n)

    for i in 1:n
        ray_mat[i+1, i] = 1
    end
    
    # Track the first row of the theoretical unimodular matrix as a separate column vector
    v1 = zero_matrix(ZZ, n_plus_1, 1)
    v1[1, 1] = 1
    
    cur_x = copy(x)
    
    # 2. Unimodular reduction (Extended Euclidean Algorithm)
    for i in 2:n_plus_1
        a = cur_x[1]
        b = cur_x[i]
        
        g, u, v = gcdx(a, b)
        
        # Apply transformation: We iterate over j, accessing contiguous memory in the column
        for j in 1:n_plus_1
            old_v1 = v1[j]
            old_vi = ray_mat[j, i-1]
            
            v1[j]           = ZZ(u * old_v1 + v * old_vi)
            ray_mat[j, i-1] = ZZ(-(b ÷ g) * old_v1 + (a ÷ g) * old_vi)
        end
        
        cur_x[1] = g
        cur_x[i] = 0
    end
    
    # 3. Incidence Matrix (See expert tip below regarding Oscar.jl)
    num_cones = 2^n_plus_1 - 1
    inc_mat = zero!(BitArray(undef, (num_cones, n_plus_1)))
    # inc_mat = Oscar.incidence_matrix(num_cones, n_plus_1)
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