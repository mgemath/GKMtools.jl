# Isotropy data for vertices and edge multiplicities for orbifold GKM graphs
# For each vertex, we store the invariants of the isotropy group and the weights of the representation on the tangent space. 
# For each edge, we store the multiplicity (the order of the isotropy group along that edge).  
struct OrbifoldVertexIsotropy <: AbstractVertexIsotropy{Vector{Int},Matrix{Int64}}
  isotropy_group::Vector{Int}  # group structure of the isotropy group at the vertex
  tangent_rep::Matrix{Int64}     # full tangent representation
end

struct OrbifoldFlagIsotropy <: AbstractFlagIsotropy{Vector{Int},Matrix{Int64},ZZMatrix}
  isotropy_group::Vector{Int}   # group structure of the isotropy group at the flag
  tangent_rep::Matrix{Int64}     # full tangent representation
  embedding::ZZMatrix        # how the isotropy group of the flag embeds into the isotropy group of the vertex
end

function order_of_isotropy_group(I::AbstractIsotropy{Vector{Int},Matrix{Int64}})
  return prod(I.isotropy_group; init=1)
end

function order_of_generic_stabilizer(V::OrbifoldVertexIsotropy, F::OrbifoldFlagIsotropy)
  return order_of_isotropy_group(V) ÷ order_of_isotropy_group(F)
end

function smooth_orbifold_vertex_isotropy_group(d)
  return OrbifoldVertexIsotropy(Int[], zeros(Int, 0, d))
end

function smooth_orbifold_flag_isotropy_group(d, k)
  return OrbifoldFlagIsotropy(Int[], zeros(Int, 0, d), identity_matrix(ZZ, k))
end

issmooth(I::AbstractIsotropy{Vector{Int},Matrix{Int64}}) = isempty(I.isotropy_group)