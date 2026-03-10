export generalized_gkm_flag

@doc raw"""
    generalized_gkm_flag(R::RootSystem, S::Vector{RootSpaceElem}) -> AbstractGKM_graph

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
    generalized_gkm_flag(R::RootSystem; indices_of_S) -> AbstractGKM_graph

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
