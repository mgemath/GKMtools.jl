
struct StackyCone <: AbstractStackyCone
  rays_matrix::ZZMatrix
end

struct WeightedProjectiveSpace <: AbstractStackyFan
  rays_matrix::ZZMatrix
  incidence_matrix::BitMatrix
  w::Vector{Int}
end

function incidence_matrix(X::AbstractStackyFan)
  return X.incidence_matrix
end

function incidence_matrix(X::AbstractStackyCone)
  return BitArray(undef, (1, size(X.rays_matrix, 1))) .|| true
end


function Oscar.cones(X::WeightedProjectiveSpace, k::Int)
  if k < 0 || k > size(X.incidence_matrix, 2)
    error("k must be between 0 and the number of rays")
  end
  indices = findall(i -> sum(X.incidence_matrix[i, :]) == k, 1:size(X.incidence_matrix, 1))
  return [cone_number(X, n) for n in indices] 
end

function Oscar.cones(X::WeightedProjectiveSpace)
  return vcat([Oscar.cones(X, k) for k in 0:n_rays(X)]...)
end

function Oscar.maximal_cones(X::WeightedProjectiveSpace)
  return Oscar.cones(X, dim(X))
end

function Oscar.maximal_cones(C::StackyCone)
  return [C]
end

function Oscar.rays(X::WeightedProjectiveSpace)
  return [X.rays_matrix[i, :] for i in 1:size(X.rays_matrix, 1)]
end

function Oscar.rays(C::StackyCone)
  return [C.rays_matrix[i, :] for i in 1:size(C.rays_matrix, 1)]
end

function Oscar.dim(X::WeightedProjectiveSpace)
  return length(X.w) - 1
end

function Oscar.dim(C::StackyCone)
  return size(C.rays_matrix, 1)
end

function Oscar.n_rays(X::WeightedProjectiveSpace)
  return size(X.incidence_matrix, 2)
end

function Oscar.n_rays(C::StackyCone)
  return size(C.rays_matrix, 1)
end

function Oscar.is_orbifold(::Union{AbstractStackyFan, AbstractStackyCone})
  return true
end

function Oscar.polarize(C::StackyCone)
  # This is a placeholder for the actual polarize function, which would compute the dual cone.
  # For simplicity, we will just return a new StackyCone with the same rays_matrix.
  polar_cone = polarize(cone(C.rays_matrix))
  return polar_cone
  # return StackyCone(C.rays_matrix)
end

function cone_number(X::WeightedProjectiveSpace, n)
  indices = findall(i -> X.incidence_matrix[n, i] == 1, 1:size(X.incidence_matrix, 2))
  return StackyCone(matrix(ZZ, [X.rays_matrix[i, :] for i in indices]))
end