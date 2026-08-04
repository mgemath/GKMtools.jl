function _check_vector_bundle_lattice(
  G::AbstractGKMGraph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
) where {R}
  @req domain(GMtoM) == lattice(G) """
  GMtoM must go from lattice(G) into M
  """

  @req codomain(GMtoM) == M """
  GMtoM must go from lattice(G) into M
  """

  return nothing
end

function _check_vector_bundle_weights(
  G::AbstractGKMGraph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  weights::Matrix{
    AbstractAlgebra.Generic.FreeModuleElem{R}
  },
) where {R}
  @req size(weights, 1) == num_vertices(G) """
  Weight matrix has wrong number of rows
  """

  for w in weights
    @req parent(w) == M "Weights need to live in M"
  end

  return nothing
end

function _trivial_fiber_representations(G::AbstractOrbifoldGKMGraph, rank_bundle::Int)
  return [
    zero_matrix(ZZ, length(G.vertex_isotropy[v].isotropy_group), rank_bundle)
    for v in 1:num_vertices(G)
  ]
end

function _check_fiber_representations(
  G::AbstractOrbifoldGKMGraph,
  fiber_reps::Vector{ZZMatrix},
  rank_bundle::Int,
)
  @req length(fiber_reps) == num_vertices(G) "fiber_reps must contain one matrix for each vertex"
  for v in 1:num_vertices(G)
    isotropy_rank = length(G.vertex_isotropy[v].isotropy_group)
    @req size(fiber_reps[v], 1) == isotropy_rank "fiber_reps[$v] has the wrong number of rows"
    @req size(fiber_reps[v], 2) == rank_bundle "fiber_reps[$v] has the wrong number of columns"
  end
  return nothing
end

@doc raw"""
    orbifold_vector_bundle(G, M, GMtoM, weights; fiber_reps=nothing)

Construct an equivariant vector bundle over an orbifold GKM graph.

The matrix `weights` stores the torus weights of the fibre over each
fixed point: `weights[v, i]` is the weight of the `i`-th line summand
over vertex `v`.

The optional `fiber_reps` stores the action of the vertex stabilizers on
the fibre.  For a rank `r` bundle, `fiber_reps[v]` is a matrix with
`length(G.vertex_isotropy[v].isotropy_group)` rows and `r` columns.  Its
`i`-th column is the character of the local isotropy group on the
`i`-th fibre summand.  If omitted, the stabilizers act trivially on the
fibres.
"""
function orbifold_vector_bundle(
  G::AbstractOrbifoldGKMGraph{R,V,F},
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}};
  fiber_reps::Union{Nothing,Vector{ZZMatrix}} = nothing,
) where {R,V,F}
  _check_vector_bundle_lattice(G, M, GMtoM)
  _check_vector_bundle_weights(G, M, weights)

  rank_bundle = size(weights, 2)
  reps = isnothing(fiber_reps) ? _trivial_fiber_representations(G, rank_bundle) : fiber_reps
  _check_fiber_representations(G, reps, rank_bundle)

  return OrbifoldGKMVectorBundle{R,V,F,typeof(G)}(G, M, GMtoM, weights, reps, build_vector_bundle_connection(G, GMtoM, weights))
end

@doc raw"""
    orbifold_line_bundle(
      G::AbstractOrbifoldGKMGraph{R,V,F},
      M::AbstractAlgebra.Generic.FreeModule{R},
      GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
      weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}};
      fiber_reps::Union{Nothing,Vector{ZZMatrix}} = nothing,
) where {R,V,F}

Construct a rank-one equivariant vector bundle over an orbifold GKM graph.
"""
function orbifold_line_bundle(
  G::AbstractOrbifoldGKMGraph{R,V,F},
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}};
  fiber_reps::Union{Nothing,Vector{ZZMatrix}} = nothing,
) where {R,V,F}
  return orbifold_vector_bundle(
    G,
    M,
    GMtoM,
    reshape(weights, length(weights), 1);
    fiber_reps = fiber_reps,
  )
end

@doc raw"""
    vector_bundle(G, M, GMtoM, weights; fiber_reps=nothing) -> AbstractGKMVectorBundle

Construct an equivariant GKM vector bundle. `G` is the base graph, `M` is
the character lattice of the fibre weights, `GMtoM` maps the character
lattice of `G` to `M`, and `weights[v, i]` is the `i`-th fibre weight at
vertex `v`.

For an orbifold base, `fiber_reps` optionally records the representations of
the vertex stabilizers. In that case this is an alias for
[`orbifold_vector_bundle`](@ref).
"""
vector_bundle(
  G::AbstractOrbifoldGKMGraph,
  M::AbstractAlgebra.Generic.FreeModule,
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism,
  weights::Matrix{<:AbstractAlgebra.Generic.FreeModuleElem};
  fiber_reps::Union{Nothing,Vector{ZZMatrix}} = nothing,
) = orbifold_vector_bundle(G, M, GMtoM, weights; fiber_reps = fiber_reps)

@doc raw"""
    line_bundle(G, M, GMtoM, weights; fiber_reps=nothing) -> AbstractGKMVectorBundle

Construct a rank-one equivariant GKM vector bundle from one fibre weight at
each vertex. For an orbifold base, `fiber_reps` optionally records the
characters of the vertex stabilizers and this function is an alias for
[`orbifold_line_bundle`](@ref).
"""
line_bundle(
  G::AbstractOrbifoldGKMGraph,
  M::AbstractAlgebra.Generic.FreeModule,
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism,
  weights::Vector{<:AbstractAlgebra.Generic.FreeModuleElem};
  fiber_reps::Union{Nothing,Vector{ZZMatrix}} = nothing,
) = orbifold_line_bundle(G, M, GMtoM, weights; fiber_reps = fiber_reps)

function Oscar.rank(
  E::AbstractGKMVectorBundle,
)::Int
  return size(E.weights, 2)
end

@doc raw"""
    baseof(V::AbstractGKMVectorBundle) -> AbstractGKM_graph

Return the base of the given GKM vector bundle.
"""
baseof(E::AbstractGKMVectorBundle) = E.base

@doc raw"""
    fiber_weight(E::AbstractGKMVectorBundle, v::Integer, i::Integer)

Return the `i`-th torus weight of the fibre of `E` over vertex `v`.
"""
function fiber_weight(
  E::AbstractGKMVectorBundle,
  v::Int,
  i::Int,
)
  return E.weights[v, i]
end

@doc raw"""
    fiber_representation(E::OrbifoldGKMVectorBundle, v::Integer) -> ZZMatrix

Return the matrix describing the action of the isotropy group at vertex `v`
on the fibre of `E`. Rows correspond to cyclic isotropy factors and columns
to the weight lines of the fibre.
"""
fiber_representation(E::OrbifoldGKMVectorBundle, v::Int) = E.fiber_reps[v]


function _tangent_bundle_weights(
  G::AbstractGKMGraph{R},
) where {R}
  tangent_rank = valency(G)

  weights = Matrix{
    AbstractAlgebra.Generic.FreeModuleElem{R}
  }(
    undef,
    num_vertices(G),
    tangent_rank,
  )

  for v in 1:num_vertices(G)
    @req length(flags(G, v)) == tangent_rank """
    Tangent bundle construction expects constant valency
    """

    for i in 1:tangent_rank
      weights[v, i] = weight(G, v, i)
    end
  end

  return weights
end

function _tangent_bundle_fiber_representations(G::AbstractOrbifoldGKMGraph)
  tangent_rank = valency(G)
  reps = Vector{ZZMatrix}(undef, num_vertices(G))
  for v in 1:num_vertices(G)
    rep = G.vertex_isotropy[v].tangent_rep
    @req size(rep, 2) == tangent_rank "Vertex tangent representation has wrong number of columns"
    reps[v] = rep
  end
  return reps
end

@doc raw"""
    tangent_bundle(G::GKMGraph)

Construct the tangent bundle of a smooth GKM graph.
"""
function Oscar.tangent_bundle(G::GKMGraph)
  M = lattice(G)
  GMtoM = hom(M, M, gens(M))
  weights = _tangent_bundle_weights(G)

  return vector_bundle(
    G,
    M,
    GMtoM,
    weights,
  )
end

@doc raw"""
    tangent_bundle(G::AbstractOrbifoldGKMGraph)

Construct the tangent bundle of an orbifold GKM graph.
"""
function Oscar.tangent_bundle(
  G::AbstractOrbifoldGKMGraph,
)
  M = lattice(G)
  GMtoM = hom(M, M, gens(M))
  weights = _tangent_bundle_weights(G)

  fiber_reps =
    _tangent_bundle_fiber_representations(G)

  return orbifold_vector_bundle(
    G,
    M,
    GMtoM,
    weights;
    fiber_reps=fiber_reps,
  )
end

function vector_bundle(
  G::GKMGraph{R,V,F},
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Matrix{
    AbstractAlgebra.Generic.FreeModuleElem{R}
  },
) where {R,V,F}
  _check_vector_bundle_lattice(G, M, GMtoM)
  _check_vector_bundle_weights(G, M, weights)

  return GKMVectorBundle{R,V,F,typeof(G)}(
    G,
    M,
    GMtoM,
    weights,
    build_vector_bundle_connection(G, GMtoM, weights),
  )
end

function line_bundle(
  G::GKMGraph{R,V,F},
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Vector{
    AbstractAlgebra.Generic.FreeModuleElem{R}
  },
) where {R,V,F}
  return vector_bundle(
    G,
    M,
    GMtoM,
    reshape(weights, length(weights), 1),
  )
end