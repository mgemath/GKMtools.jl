export generalized_gkm_flag

@doc raw"""
    generalized_gkm_flag(R::RootSystem, S::Vector{RootSpaceElem}) -> GKMGraph

Given a root system ``R`` and a subset ``S`` of the set of simple roots, it constructs the 
GKM graph of the generalized flag variety ``G/P``. Here ``G`` is the simply-connected complex Lie group 
with root system ``R``, and ``P`` is the parabolic subgroup with root system ``S``.
If ``S`` is empty, it construct ``G/B`` where ``B`` is a Borel subgroup.
The vertices of ``G/P`` correspond to the cosets ``W/W_P`` where ``W`` (resp., ``W_P``) is the Weyl group
of ``G`` (resp., ``P``). The label of a vertex is the unique element of minimal length in the corresponding coset.

!!! note
    The character group is of type free ``\mathbb{Z}``-module if ``R`` is of type ``A, B, C, D, G``.
    It is a free ``\mathbb{Q}``-module if ``R`` is of type ``E`` or ``F``.

!!! warning
    Computing this function with root systems of very large Weyl groups may be slow.

# Examples
```jldoctest
julia> A1xA1 = root_system([(:A, 1), (:A, 1)])
Root system of rank 2
  of type A1 x A1

julia> generalized_gkm_flag(A1xA1)
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (-1, 1, 0, 0)
s2 -> id => (0, 0, -1, 1)
s1*s2 -> s1 => (0, 0, -1, 1)
s1*s2 -> s2 => (-1, 1, 0, 0)

julia> RC3 = root_system(:C, 3)
Root system of rank 3
  of type C3

julia> gp1 = generalized_gkm_flag(RC3);

julia> rank_torus(gp1)
3

julia> R = root_system([(:A, 1), (:G, 2)])
Root system of rank 3
  of type A1 x G2

julia> S = [simple_roots(R)[3]];

julia> gp2 = generalized_gkm_flag(R, S);

julia> rank_torus(gp2)
5

```
"""
function generalized_gkm_flag(R::RootSystem, S::Vector{RootSpaceElem})
  @req all(sr -> sr in simple_roots(R), S) "S must be a set of simple roots of R"

  return generalized_gkm_flag(R, findall(j -> simple_root(R, j) in S, 1:rank(R)))
end

@doc raw"""
    generalized_gkm_flag(R::RootSystem; indices_of_S) -> GKMGraph

Same as before, but indicating the indices of the roots in ``S`` instead of the roots itself.

# Examples
```jldoctest
julia> R = root_system(matrix(ZZ, [2 -1 -2; -1 2 0; -1 0 2]))
Root system of rank 3
  of type C3 (with non-canonical ordering of simple roots)

julia> gp1 = generalized_gkm_flag(R, 2:3);

julia> valency(gp1)
7

julia> gp2 = generalized_gkm_flag(R, [1,2]);

julia> rank_torus(gp2)
3

```
"""
function generalized_gkm_flag(R::RootSystem, indices_of_S::AbstractVector{<:Integer}=Int[])
  _check_consistency(R, indices_of_S)

  # 1. Create WP
  Weyl = weyl_group(R)
  WP = _WP(R, indices_of_S)

  # 2. Geometry (Cosets & Graph)
  cosets, reprs = _cosets_and_repr(Weyl, WP)
  gen_matrix, type_of_graph = _gen_matrix_and_type_of_graph(R)

  # 3. Construct and return GP
  GP = _generalized_gkm_flag(R, cosets, reprs, WP, gen_matrix, type_of_graph)

  return GP
end

function _dimension_ambient(RT::Tuple{Symbol,Int64})::Int64
  if RT[1] in (:A, :G)
    return RT[2] + 1
  elseif RT[1] == :E
    return 8
  end

  return RT[2]
end

function _generator_matrix(RT::Tuple{Symbol,Int64})::QQMatrix
  n_rows = RT[2]
  n_cols = _dimension_ambient(RT)

  if RT[1] == :E # following Hum75 convention, pag 64
    M = zero_matrix(QQ, n_rows, n_cols)
    foreach(i -> M[1, i] = (i == 1 || i == 8) ? QQ(1)//QQ(2) : QQ(-1)//QQ(2), 1:n_cols)
    M[2, 1] = QQ(1)
    M[2, 2] = QQ(1)
    for i in 3:n_rows
      M[i, i - 2] = QQ(-1)
      M[i, i - 1] = QQ(1)
    end
    return M
  end

  M = diagonal_matrix(QQ(1), n_rows, n_cols)

  if RT[1] == :G  # following Hum75 convention
    M[1, 2] = QQ(-1)
    M[2, 1] = QQ(-2)
    M[2, 3] = QQ(1)
    return M
  end

  for i in 2:n_cols
    M[i - 1, i] = QQ(-1)
  end

  if RT[1] in (:A, :B) # # following Hum75 & Wikipedia convention
  else
    if RT[1] == :C # following Hum75 & Wikipedia convention
      M[n_rows, n_cols] = QQ(2)
    elseif RT[1] == :D # following Hum75 & Wikipedia convention
      M[n_rows, n_cols - 1] = QQ(1)
    elseif RT[1] == :F # following Wikipedia convention
      M[3, 4] = QQ(0)
      foreach(i -> M[4, i] = QQ(-1)//QQ(2), 1:4)
    end
  end

  return M
end

function _connection_from_cartan_data(core::GKMCombinatorialData{R}, cartan_data) where {R}
  bundle_transport = Dict{Edge,Vector{Int}}()
  bundle_coefficients = Dict{Edge,Vector{R}}()
  flag_edges = [Vector{Edge}(undef, length(flags(core, v))) for v in vertices(core)]
  for edge in edges(graph(core))
    source_flag, target_flag = core.edge_flags[edge]
    flag_edges[src(edge)][source_flag] = edge
    flag_edges[dst(edge)][target_flag] = reverse(edge)
  end

  for base_edge in edges(graph(core)), e in (base_edge, reverse(base_edge))
    source_weights = flags(core, src(e))
    target_weights = flags(core, dst(e))
    target_lookup = Dict(flag.weight => i for (i, flag) in enumerate(target_weights))
    image = Vector{Int}(undef, length(source_weights))
    coeffs = Vector{R}(undef, length(source_weights))
    edge_weight = weight(core, e)

    for i in eachindex(source_weights)
      source_flag_edge = flag_edges[src(e)][i]
      coefficient = R(cartan_data[(e, source_flag_edge)])
      target_weight = source_weights[i].weight - coefficient * edge_weight
      j = get(target_lookup, target_weight, 0)
      j == 0 && throw(ArgumentError("in-place connection has no target for flag $i along $e"))
      image[i], coeffs[i] = j, coefficient
    end

    bundle_transport[e] = image
    bundle_coefficients[e] = coeffs
  end

  return Connection{R}(bundle_transport, bundle_coefficients, "Lie")
end

function _generalized_gkm_flag(R, cosets, reprs, WP, gen_matrix, ::Type{C}) where {C}
  coset_map = Dict(element => index for (index, coset) in enumerate(cosets) for element in coset)
  WP_set = Set(WP)
  g = Graph{Undirected}(length(reprs))
  M = free_module(parent(zero(C)), ncols(gen_matrix))
  basis = gens(M)
  flags = [FlagWeight{C}[] for _ in reprs]
  edge_flags = Dict{Edge,Tuple{Int,Int}}()
  roots = Dict{Edge,RootSpaceElem}()
  positive = [root for root in positive_roots(R) if reflection(root) ∉ WP_set]

  for (i, omega) in enumerate(reprs)
    inverse_omega = inv(omega)
    for root in positive
      j = coset_map[omega * reflection(root)]
      roots[Edge(i, j)] = root
      j > i || continue

      coeffs = matrix(parent(zero(C)), Oscar.coefficients(root * inverse_omega) * gen_matrix)
      axial_weight = -sum(k -> coeffs[k] * basis[k], eachindex(basis); init=zero(M))
      add_edge!(g, j, i)
      e = Edge(j, i)
      push!(flags[j], FlagWeight{C}(axial_weight))
      push!(flags[i], FlagWeight{C}(-axial_weight))
      edge_flags[e] = (length(flags[j]), length(flags[i]))
    end
  end

  labels = GeneralizedFlagVertex[
    GeneralizedFlagVertex(replace(repr(representative), " " => ""), representative)
    for representative in reprs
  ]
  core = GKMCombinatorialData{C,GeneralizedFlagVertex,FlagWeight{C}}(
    g, M, labels, flags, edge_flags,
  )

  cartan_data = Dict{Tuple{Edge,Edge},ZZRingElem}()
  for v in vertices(g)
    neighbors = collect(all_neighbors(g, v))
    for w in neighbors
      alpha = roots[Edge(v, w)]
      for u in neighbors
        beta = roots[Edge(v, u)]
        cartan_data[(Edge(v, w), Edge(v, u))] =
          ZZ(2 * dot(beta, alpha) // dot(alpha, alpha))
      end
    end
  end

  connection = _connection_from_cartan_data(core, cartan_data)
  cohomology = create_cohomology(rank(M), length(reprs))
  return GKMGraph{C,GeneralizedFlagVertex,FlagWeight{C}}(
    core, connection, cohomology, nothing, nothing,
  )
end
function _WP(R, indices_of_S)

  if isempty(indices_of_S)
    return [one(weyl_group(R))] # WP is the whole Weyl group if S is empty
  end

  # No changes to logic, just cleaner iteration
  levi_root_system = _levi_subroot_system(R, indices_of_S)

  return [
    prod(i -> reflection(simple_root(R, indices_of_S[i])), word(a); init=one(weyl_group(R))) for
    a in weyl_group(levi_root_system)
  ]
end

function _cosets_and_repr(Weyl, WP)
  # Optimized lookup using Set, but strictly preserving the Left Coset logic (b * WP)

  total_cosets = div(order(Weyl), length(WP))
  cosets = Vector{Vector{eltype(Weyl)}}(undef, Int(total_cosets))
  reprs = Vector{eltype(Weyl)}(undef, Int(total_cosets))

  # WP is the first coset (Identity * WP)
  cosets[1] = collect(WP)
  reprs[1] = one(Weyl) # Identity

  # Fast lookup to check if element is already covered
  covered = Set{eltype(Weyl)}(WP)

  index = 1

  for b in Weyl
    if b in covered
      continue
    end

    index += 1
    # Original logic: cosets[index] = b .* cosets[index] (where cosets[index] starts as WP)
    # We replicate this: new_coset = b * WP
    new_coset = [b * w for w in WP]

    cosets[index] = new_coset

    # Find representative (shortest length)
    reprs[index] = reduce((x, y) -> length(x) <= length(y) ? x : y, new_coset)

    # Mark these as covered
    union!(covered, new_coset)

    if index == total_cosets
      break
    end
  end

  return cosets, reprs
end

function _gen_matrix_and_type_of_graph(R::RootSystem)
  (fams, ordering) = root_system_type_with_ordering(R)

  type_of_graph = any(fam -> fam[1] in (:E, :F), fams) ? QQFieldElem : ZZRingElem

  return AbstractAlgebra.perm(ordering) *
         block_diagonal_matrix([_generator_matrix(fam) for fam in fams]),
  type_of_graph
end

