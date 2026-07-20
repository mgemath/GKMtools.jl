
"""A simplicial stacky cone whose rows are the stacky ray generators."""
struct StackyCone <: AbstractStackyCone
  rays_matrix::ZZMatrix
  primitive_rays_matrix::ZZMatrix
  ambient_rank::Int
  ray_indices::Vector{Int}
end

function _primitive_rays(rays_matrix::ZZMatrix)
  primitive = zero_matrix(ZZ, nrows(rays_matrix), ncols(rays_matrix))
  for i in 1:nrows(rays_matrix)
    entries = [Int(rays_matrix[i, j]) for j in 1:ncols(rays_matrix)]
    divisor = foldl(gcd, abs.(entries); init=0)
    divisor == 0 && throw(ArgumentError("stacky ray $i must be nonzero"))
    for j in 1:ncols(rays_matrix)
      primitive[i, j] = div(rays_matrix[i, j], divisor)
    end
  end
  return primitive
end

function StackyCone(rays_matrix::ZZMatrix; ambient_rank::Integer=Int(ncols(rays_matrix)), ray_indices::AbstractVector{<:Integer}=collect(1:Int(nrows(rays_matrix))))
  nr, ar = Int(nrows(rays_matrix)), Int(ambient_rank)
  ar == ncols(rays_matrix) || throw(ArgumentError("ambient_rank must equal the number of matrix columns"))
  length(ray_indices) == nr || throw(ArgumentError("one ray index is required for each stacky generator"))
  length(unique(ray_indices)) == nr || throw(ArgumentError("ray indices must be distinct"))
  primitive = _primitive_rays(rays_matrix)
  length(unique([Tuple(primitive[i, :]) for i in 1:nr])) == nr || throw(ArgumentError("a cone cannot contain the same geometric ray twice"))
  rank(rays_matrix) == nr || throw(ArgumentError("an orbifold stacky cone must be simplicial"))
  return StackyCone(rays_matrix, primitive, ar, Int.(ray_indices))
end

StackyCone(rays::AbstractMatrix{<:Integer}; kwargs...) = StackyCone(matrix(ZZ, rays); kwargs...)

"""A validated orbifold stacky fan, with maximal cones given by ray indices."""
struct StackyFan <: AbstractStackyFan
  rays_matrix::ZZMatrix
  primitive_rays_matrix::ZZMatrix
  maximal_cone_indices::Vector{Vector{Int}}
  ambient_rank::Int
  generic_stabilizer::Vector{Int}
end

function StackyFan(rays_matrix::ZZMatrix, maximal_cones::AbstractVector{<:AbstractVector{<:Integer}}; generic_stabilizer::AbstractVector{<:Integer}=Int[])
  nr, ar = Int(nrows(rays_matrix)), Int(ncols(rays_matrix))
  nr > 0 || throw(ArgumentError("a stacky fan needs at least one ray"))
  cones = [Int.(c) for c in maximal_cones]
  isempty(cones) && throw(ArgumentError("a stacky fan needs at least one maximal cone"))
  all(>(1), generic_stabilizer) || throw(ArgumentError("generic stabilizer invariants must be greater than one"))
  primitive = _primitive_rays(rays_matrix)
  length(unique([Tuple(primitive[i, :]) for i in 1:nr])) == nr || throw(ArgumentError("stacky generators must lie on distinct geometric rays"))
  for (i, indices) in enumerate(cones)
    length(indices) == ar || throw(ArgumentError("maximal cone $i must have $ar rays"))
    length(unique(indices)) == length(indices) || throw(ArgumentError("maximal cone $i contains a ray more than once"))
    all(j -> 1 <= j <= nr, indices) || throw(ArgumentError("maximal cone $i contains an invalid ray index"))
    rank(rays_matrix[indices, :]) == ar || throw(ArgumentError("maximal cone $i is not full-dimensional and simplicial"))
  end
  return StackyFan(rays_matrix, primitive, cones, ar, Int.(generic_stabilizer))
end

StackyFan(rays::AbstractMatrix{<:Integer}, maximal_cones; kwargs...) = StackyFan(matrix(ZZ, rays), maximal_cones; kwargs...)

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

function incidence_matrix(X::StackyFan)
  result = falses(length(X.maximal_cone_indices), n_rays(X))
  for (i, indices) in enumerate(X.maximal_cone_indices)
    result[i, indices] .= true
  end
  return result
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

function Oscar.maximal_cones(X::StackyFan)
  return [StackyCone(X.rays_matrix[indices, :]; ambient_rank=X.ambient_rank, ray_indices=indices) for indices in X.maximal_cone_indices]
end

function Oscar.maximal_cones(C::StackyCone)
  return [C]
end

function Oscar.rays(X::WeightedProjectiveSpace)
  return [X.rays_matrix[i, :] for i in 1:size(X.rays_matrix, 1)]
end

Oscar.rays(X::StackyFan) = [X.rays_matrix[i, :] for i in 1:nrows(X.rays_matrix)]

function Oscar.rays(C::StackyCone)
  return [C.rays_matrix[i, :] for i in 1:size(C.rays_matrix, 1)]
end

function Oscar.dim(X::WeightedProjectiveSpace)
  return length(X.w) - 1
end

function Oscar.dim(C::StackyCone)
  return Int(rank(C.rays_matrix))
end

Oscar.dim(X::StackyFan) = X.ambient_rank

function Oscar.n_rays(X::WeightedProjectiveSpace)
  return size(X.incidence_matrix, 2)
end

function Oscar.n_rays(C::StackyCone)
  return size(C.rays_matrix, 1)
end

Oscar.n_rays(X::StackyFan) = Int(nrows(X.rays_matrix))

Oscar.is_orbifold(C::StackyCone) = rank(C.rays_matrix) == n_rays(C)
Oscar.is_orbifold(X::StackyFan) = all(C -> dim(C) == dim(X) && n_rays(C) == dim(X), maximal_cones(X))
Oscar.is_orbifold(X::WeightedProjectiveSpace) = all(C -> dim(C) == dim(X) && n_rays(C) == dim(X), maximal_cones(X))

function Oscar.polarize(C::StackyCone)
  # This is a placeholder for the actual polarize function, which would compute the dual cone.
  # For simplicity, we will just return a new StackyCone with the same rays_matrix.
  polar_cone = polarize(cone(C.rays_matrix))
  return polar_cone
  # return StackyCone(C.rays_matrix)
end

function cone_number(X::WeightedProjectiveSpace, n)
  indices = findall(i -> X.incidence_matrix[n, i] == 1, 1:size(X.incidence_matrix, 2))
  return StackyCone(matrix(ZZ, [X.rays_matrix[i, :] for i in indices]); ambient_rank=dim(X), ray_indices=indices)
end

function Base.show(io::IO, X::WeightedProjectiveSpace)
  println(io, "Weighted Projective Space P(", join(X.w, ", "), ")")
end