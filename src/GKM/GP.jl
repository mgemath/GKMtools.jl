@doc raw"""
    generalized_gkm_flag(R::RootSystem, S::Vector{RootSpaceElem}) -> GKM_graph

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

  s = simple_roots(R)

  return _generalized_gkm_flag(R, findall(j -> s[j] in S, eachindex(s)))
end

@doc raw"""
    generalized_gkm_flag(R::RootSystem; indices_of_S) -> GKM_graph

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
function generalized_gkm_flag(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}=Int64[])

  @req isempty(indices_of_S) || (minimum(indices_of_S) > 0 && maximum(indices_of_S) <= rank(R)) "indices of S out of range"

  return _generalized_gkm_flag(R, collect(indices_of_S))
end

function _generalized_gkm_flag(R::RootSystem, indices_of_S::Vector{Int64})

  ## Create WP
  Weyl = weyl_group(R) # Weyl group of the root system
  WP = [one(Weyl)] # embedding of W_P into W, called WP

  if !isempty(indices_of_S)

    #Cartan matrix of P
    cartan_submatrix = cartan_matrix(R)[indices_of_S, indices_of_S]

    #embedding of W_P into W, called WP
    s = simple_roots(R)
    WP = [prod(i -> reflection(s[indices_of_S[i]]), word(a); init = one(Weyl)) for a in weyl_group(cartan_submatrix)]
  end
  
  (fams, ordering) = root_system_type_with_ordering(R)

  type_of_graph = any(fam -> fam[1] in (:E, :F), fams) ? QQFieldElem : ZZRingElem
  
  cosets = [WP for _ in 1:div(order(Weyl), length(WP))]
  reprs = [one(Weyl) for _ in 1:length(cosets)]
  index = 1
  for b in Weyl
    any(i -> b in cosets[i], 1:index) && continue
    index += 1
    cosets[index] = b .* cosets[index]
    # cosets[index] = [a * b for a in WP]
    reprs[index] = reduce((x, y) -> length(x) <= length(y) ? x : y,  cosets[index])
    index == length(cosets) && break
  end

  #println("cosets: $cosets")
  #println("reprs: $reprs")

  gen_matrix = AbstractAlgebra.perm(ordering)*block_diagonal_matrix([_generator_matrix(fam) for fam in fams])
  labs = [replace(repr(r), " " => "") for r in reprs]# repr.(reprs)
  # println("gen_matrix = $gen_matrix")
  # labs = ["$i" for i in 1:length(reprs)]# repr.(reprs)
  g = Graph{Undirected}(length(labs))
  M = free_module(parent(zero(type_of_graph)), n_columns(gen_matrix))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{type_of_graph}}()

  # println("cosets[1] = $(cosets[1])")
  # println("cosets[2] = $(cosets[2])")
  # println("cosets[3] = $(cosets[3])")
  get_root = Dict{Edge, RootSpaceElem}()

  for i in 1:length(reprs)
    omega = reprs[i]
    for t in positive_roots(R)  #Delta_G

      # Remove roots from S
      reflection(t) in WP && continue  #now I run on Delta_G,K
      
      # A representative of the coset which we are connecting to
      # new_rep = reflection(t)*r1
      
      # new_rep = reflection(t*r1)
      # new_rep = reflection(t)*r1
      new_rep = omega * reflection(t)   # w*sigma_t
      # new_rep = reflection(t) * omega
      # j = findfirst(index -> (index > i) && (new_rep in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 
      # j = findfirst(index ->  (index != i && new_rep in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 

      # Find the representative of the new coset smallest length
      j = findfirst(var ->  var != i && new_rep in cosets[var], 1:length(reprs)) #change index != i if you want double ways 
      # j === nothing && error() # a representative must existk
      j === nothing && continue

      # Let us deal with the weights, we add the weight to W only if (j > i)

      
      if j > i
        # We orient the edge from maximum(i,j) to minimum(i,j)
        add_edge!(g, j, i)

        sign = -1 # This is because the calculated weight is for Edge(i, j), but then we set it for Edge(j, i).

        # The weight:-----------------------
        # Explanation for inv(omega):
        # OSCAR only supports right mutiplication by Weyl group elements.
        # to right multiply a root by w = s1 s2 ... sn, it reflects the root first by s1, then by s2, ..., then by sn.
        # Left multiplication by w should be reflecting the root by sn first, then by s(n-1), ..., then by s2, then by s1.
        # Thus, omega * t = t * inv(omega).
        vec = matrix(parent(zero(type_of_graph)), coefficients(t*inv(omega))*gen_matrix)    #t*r1 = w(root)  
        W[Edge(j, i)] = sign*sum(_i -> vec[_i]*gens(M)[_i], 1:rank(M))
        
        #println("t=$t, omega=$omega, new_rep = $new_rep, sign  = $sign")
        #println("t*inv(omega)=$(t*inv(omega))")
        #println("vec = $vec")
        #println("W[Edge($j, $i)]=$(W[Edge(j, i)])")
        #println("")
        
      end

      # if i==3
      #   println("i=$i, j=$j, t*omega=$(t*omega), vec=$(matrix(parent(zero(type_of_graph)), coefficients(t*omega)*gen_matrix))")
      # end

      # if j==3
      #   println("i=$i, t*omega=$(t*omega), vec=$(matrix(parent(zero(type_of_graph)), coefficients(t*omega)*gen_matrix))")
      # end

      # get_root[Edge(i, j)] = is_positive_root(t*r1) ? t*r1 : -t*r1
      get_root[Edge(i, j)] = t
      # println(Edge(i, j), " ", t)
      # println("r1=$r1, reprs[j]=$(reprs[j]), new_rep=$new_rep, t=$t, t*r1=$(t*r1), t*r1^-1=$(t*((r1)^(-1)))")
    end
  end
  #println(get_root)

  ## construct connection
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} = Dict{Tuple{Edge, Edge}, ZZRingElem}()
  # a = Dict()
  for _v in vertices(g)
    for _w in all_neighbors(g, _v)

      alpha = get_root[Edge(_v, _w)]


      for _u in all_neighbors(g, _v)
        # e = (_v, _w)
        # e'= (_v, _u)

        beta = get_root[Edge(_v, _u)]
        
        a[(Edge(_v, _w), Edge(_v, _u))] = ZZ(2*dot(beta, alpha)//dot(alpha, alpha))
      end
    end
  end

  # for e in edges(g)

  #   alpha = get_root[e]
  #   E = (src(e), dst(e))

  #   for i in 1:2
  #     _v = E[i]
  #     for _w in all_neighbors(g, _v)

  #       ei = Edge(_v, _w)
  #       # beta = get_root[ei in edges(g) ? ei : reverse(ei)]
  #       beta = get_root[ei]
        
  #       # Daniel: I removed the minus sign here, as otherwise we got an error
  #       # for A1. (We always need to have a[(e,e)] = 2 and it was -2.).
  #       # However, we still get an error for full flags in type A2!
  #       # So there is some inconsistency with the weights and the a's.
  #       # a[(Edge(E[i], E[3-i]), ei)] = ZZ(2*dot(beta, alpha)//dot(alpha, alpha))
  #       a[(Edge(E[i], E[3-i]), ei), alpha, beta] = ZZ(2*dot(beta, alpha)//dot(alpha, alpha))
  #     end
  #   end
  # end
  GP = gkm_graph(g, labs, M, W)
  
  #println(a)
  con = build_gkm_connection(GP, a)
  #println(con); 
  set_connection!(GP, con)
  return GP
end

function _dimension_ambient(RT::Tuple{Symbol, Int64})::Int64

  if RT[1] in (:A, :G)
    return RT[2] + 1
  elseif RT[1] == :E
    return 8
  end

  return RT[2] 

end

function _generator_matrix(RT::Tuple{Symbol, Int64})::QQMatrix

  n_rows = RT[2]
  n_cols = _dimension_ambient(RT)

  if RT[1] == :E # following Hum75 convention, pag 64
    M = zero_matrix(QQ, n_rows, n_cols)
    foreach(i -> M[1, i] = (i==1 || i==8) ? QQ(1)//QQ(2) : QQ(-1)//QQ(2), 1:n_cols)
    M[2, 1] = QQ(1)
    M[2, 2] = QQ(1)
    for i in 3:n_rows
      M[i, i-2] = QQ(-1)
      M[i, i-1] = QQ(1)
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
    M[i-1, i] = QQ(-1)
  end

  if RT[1] in (:A, :B) # # following Hum75 & Wikipedia convention
  else
    if RT[1] == :C # following Hum75 & Wikipedia convention
      M[n_rows, n_cols] = QQ(2)
    elseif RT[1] == :D # following Hum75 & Wikipedia convention
      M[n_rows, n_cols-1] = QQ(1)
    elseif RT[1] == :F # following Wikipedia convention
      M[3, 4] = QQ(0)
      foreach(i -> M[4, i] = QQ(-1)//QQ(2), 1:4)
    end
  end

  return M
end
