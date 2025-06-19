export get_bruhat_order_of_generalized_flag, generalized_gkm_schubert, isrationally_smooth

struct BruhatOrder
  order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}
  reversed::Bool # true for descending order
  labels::Vector{String}
end

function _add_to_order(G::AbstractGKM_graph, _order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}, elem::Int64, increment::Int64, addall::Bool = false)

    ## Create subgraph for elem
    elem_in_dict = (elem, count('s', G.labels[elem]))
    _order[elem_in_dict] = Tuple{Int64, Int64}[]
    
    
    for v in neighbors(G.g, elem)
      if (count('s', G.labels[v]) == count('s', G.labels[elem]) + increment) || (addall && count('s', G.labels[v]) < count('s', G.labels[elem]))
        push!(_order[elem_in_dict], (v, count('s', G.labels[v])))
        _add_to_order(G, _order, v, increment, addall)
      end
    end
  
  end

function _create_order(G::AbstractGKM_graph, starting_elem::Int64, rev::Bool, addall::Bool = false)::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}
  
  _order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}} = Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}()

  _add_to_order(G, _order, starting_elem, (rev ? -1 : 1), addall)

  return _order
end

function _get_bruhat_order_of_generalized_flag(G::AbstractGKM_graph, rev::Bool = false)::BruhatOrder

  starting_elem::Int64 = rev ? reduce((x,y) -> count('s', G.labels[x]) >= count('s', G.labels[y]) ? x : y, 1:n_vertices(G.g)) : 1

  _order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}} = _create_order(G, starting_elem, rev)

  return BruhatOrder(_order, rev, G.labels)

end


@doc raw"""
    get_bruhat_order_of_generalized_flag(R::RootSystem, S::Vector{RootSpaceElem}; descending::Bool=true) -> BruhatOrder

It returns the Bruhat order of the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtest.generalized_gkm_flag).
If `descending` is `true`, the Bruhat order is given from the elements of maximal length to the smallest.

# Examples
```jldoctest
julia> R = root_system(:A, 3);

julia> S = simple_roots(R);

julia> get_bruhat_order_of_generalized_flag(R, S[1:2])
Bruhat order in descending order:
Length: 3
  s1*s2*s3 => ["s2*s3"]
Length: 2
  s2*s3 => ["s3"]
Length: 1
  s3 => ["id"]
Length: 0
  id

```
"""
function get_bruhat_order_of_generalized_flag(R::RootSystem, S::Vector{RootSpaceElem}; descending::Bool=true)::BruhatOrder

  base = generalized_gkm_flag(R, S)
  
  return _get_bruhat_order_of_generalized_flag(base, descending)
end


@doc raw"""
    get_bruhat_order_of_generalized_flag(R::RootSystem, indices_of_S; descending::Bool=true) -> BruhatOrder

Same as before, but indicating the indices of the roots in ``S`` instead of the roots itself.

# Examples
```jldoctest
julia> R = root_system(:A, 3);

julia> get_bruhat_order_of_generalized_flag(R, [1]; descending = false)
Bruhat order in ascending order:
Length: 0
  id => ["s2", "s3"]
Length: 1
  s2 => ["s1*s2", "s3*s2", "s2*s3"]
  s3 => ["s3*s2", "s2*s3"]
Length: 2
  s2*s3 => ["s2*s3*s2", "s1*s2*s3"]
  s1*s2 => ["s1*s3*s2", "s1*s2*s3"]
  s3*s2 => ["s1*s3*s2", "s2*s3*s2"]
Length: 3
  s1*s2*s3 => ["s1*s2*s3*s2"]
  s2*s3*s2 => ["s2*s1*s3*s2", "s1*s2*s3*s2"]
  s1*s3*s2 => ["s2*s1*s3*s2", "s1*s2*s3*s2"]
Length: 4
  s2*s1*s3*s2 => ["s1*s2*s1*s3*s2"]
  s1*s2*s3*s2 => ["s1*s2*s1*s3*s2"]
Length: 5
  s1*s2*s1*s3*s2

```
"""
function get_bruhat_order_of_generalized_flag(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}=Int64[]; descending::Bool=true)::BruhatOrder

  base = generalized_gkm_flag(R, indices_of_S)
    
  return _get_bruhat_order_of_generalized_flag(base, descending)
end


function Base.show(io::IO, BO::BruhatOrder)

  if Oscar.is_terse(io)
      # no nested printing
    print(io, "Bruhat order")
  else
      # nested printing allowed, preferably terse
      print(io, "Bruhat order in "*(BO.reversed ? "descending" : "ascending" )*" order")
  end
end
  
  # detailed show
function Base.show(io::IO, ::MIME"text/plain", BO::BruhatOrder)
  
  print(io, "Bruhat order in "*(BO.reversed ? "descending" : "ascending" )*" order:")

  max_length = maximum(lab -> lab[2], keys(BO.order))

  seq = BO.reversed ? (max_length:-1:0) : (0:max_length)
  
  for ord in seq
    print(io, "\nLength: $ord")
    for lab in keys(BO.order)
      if lab[2] == ord
        if ord == last(seq)
          print(io, "\n  $(BO.labels[lab[1]])")
        else
          print(io, "\n  $(BO.labels[lab[1]]) => $([BO.labels[a[1]] for a in BO.order[lab]])")
        end
      end
    end
  end
  
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, indices_of_S::Vector{RootSpaceElem}, pt::String) -> AbstractGKM_subgraph

Let ``G`` be the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtest.generalized_gkm_flag).
This functions returns the subgraph of the variety ``G`` given by all Schubert cells corresponding to the points less or equal to `pt` in the Bruhat order.

# Examples
```jldoctest generalized_gkm_schubert
julia> R = root_system(:B, 2);

julia> generalized_gkm_schubert(R, "s1*s2")
GKM subgraph of:
GKM graph with 8 nodes, valency 4 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> s1 => (0, -1)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s2*s1 => (-1, 1)
s2 -> id => (0, -1)
s2 -> s2*s1 => (1, 1)
s1*s2 -> s1 => (-1, 0)
s1*s2 -> s1*s2*s1 => (1, 1)
s1*s2 -> s2 => (-1, 1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2*s1 => (-1, 0)
s2*s1*s2 -> s1*s2 => (0, -1)
s1*s2*s1*s2 -> s1 => (-1, -1)
s1*s2*s1*s2 -> s1*s2*s1 => (0, -1)
s1*s2*s1*s2 -> s2 => (-1, 0)
s1*s2*s1*s2 -> s2*s1*s2 => (-1, 1)
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (-1, 1)
s2 -> id => (0, -1)
s1*s2 -> s1 => (-1, 0)
s1*s2 -> s2 => (-1, 1)
```

As before, the subset S can be a subset of simple roots or a subset of indices.
```jldoctest generalized_gkm_schubert
julia> generalized_gkm_schubert(R, [1], "s2")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s2 -> id => (0, -1)
s1*s2 -> id => (-1, 0)
s1*s2 -> s2 => (-1, 1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2 => (-1, 0)
s2*s1*s2 -> s1*s2 => (0, -1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s2 -> id => (0, -1)

julia> S = simple_roots(R);

julia> generalized_gkm_schubert(R, [S[2]], "s1")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> id => (-1, -1)
s2*s1 -> s1 => (0, -1)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s1 => (-1, -1)
s1*s2*s1 -> s2*s1 => (-1, 1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

```
"""
function generalized_gkm_schubert(R::RootSystem, S::Vector{RootSpaceElem}, pt::String)::AbstractGKM_subgraph

  base = generalized_gkm_flag(R, S)
  @req pt in base.labels "$pt is not a valid label"

  return _generalized_gkm_schubert(base, pt, true)
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}, pt::String)  -> AbstractGKM_subgraph
"""
function generalized_gkm_schubert(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}, pt::String)::AbstractGKM_subgraph

  base = generalized_gkm_flag(R, indices_of_S)
  @req pt in base.labels "$pt is not a valid label"
  
  return _generalized_gkm_schubert(base, pt, true)
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, pt::String) -> AbstractGKM_subgraph
"""
function generalized_gkm_schubert(R::RootSystem, pt::String)::AbstractGKM_subgraph

    base = generalized_gkm_flag(R, Int64[])
    @req pt in base.labels "$pt is not a valid label"
    
    return _generalized_gkm_schubert(base, pt, true)
  end

function _generalized_gkm_schubert(base::AbstractGKM_graph, pt::String, rev::Bool=true)::AbstractGKM_subgraph

  starting_elem = findfirst(v -> base.labels[v] == pt, 1:n_vertices(base.g))
  suborder = _create_order(base, starting_elem, rev)
  

  vertices::Vector{Int64} = Int64[starting_elem]
  for lab in keys(suborder)
    vertices = vcat(vertices, [a[1] for a in suborder[lab]])
  end

  return gkm_subgraph_from_vertices(base, vertices)
end

@doc raw"""
    schubert_class(schubert::AbstractGKM_subgraph, BO::BruhatOrder, label::String)

Return the Poincare dual (as subvariety of the Schubert variety `schubert`) of the Schubert variety given by `label`.

# Arguments
- `schubert::AbstractGKM_subgraph`: The Schubert variety given as GKM subgraph object of the corresponding generalized partial flag variety (as returned by `generalized_gkm_schubert`).
- `BO::BruhatOrder`: The Bruhat order on the generalized partial flag variety (as returned by `get_bruhat_order_of_generalized_flag`).
- `label:String`: The label of the Weyl group element or coset defining the Schubert variety as subvariety of the generalized partial flag variety.

# Example
```jldoctest schubert_class
julia> R = root_system(:A, 3);

julia> S = simple_roots(R);

julia> BO = get_bruhat_order_of_generalized_flag(R, S[1:1]);

julia> schubert = generalized_gkm_schubert(R, S[1:1], "s1*s2*s3")
GKM subgraph of:
GKM graph with 12 nodes, valency 5 and axial function:
s2 -> id => (0, -1, 1, 0)
s1*s2 -> id => (-1, 0, 1, 0)
s1*s2 -> s2 => (-1, 1, 0, 0)
s3*s2 -> s2 => (0, 0, -1, 1)
s1*s3*s2 -> s1*s2 => (0, 0, -1, 1)
s1*s3*s2 -> s3*s2 => (-1, 1, 0, 0)
s2*s1*s3*s2 -> s1*s2 => (0, -1, 0, 1)
s2*s1*s3*s2 -> s1*s3*s2 => (0, -1, 1, 0)
s1*s2*s1*s3*s2 -> s2 => (-1, 0, 0, 1)
s1*s2*s1*s3*s2 -> s3*s2 => (-1, 0, 1, 0)
s1*s2*s1*s3*s2 -> s2*s1*s3*s2 => (-1, 1, 0, 0)
s2*s3*s2 -> id => (0, -1, 0, 1)
s2*s3*s2 -> s3*s2 => (0, -1, 1, 0)
s2*s3*s2 -> s2*s1*s3*s2 => (1, 0, -1, 0)
s1*s2*s3*s2 -> id => (-1, 0, 0, 1)
s1*s2*s3*s2 -> s1*s3*s2 => (-1, 0, 1, 0)
s1*s2*s3*s2 -> s1*s2*s1*s3*s2 => (0, 1, -1, 0)
s1*s2*s3*s2 -> s2*s3*s2 => (-1, 1, 0, 0)
s3 -> id => (0, 0, -1, 1)
s3 -> s3*s2 => (0, 1, 0, -1)
s3 -> s1*s3*s2 => (1, 0, 0, -1)
s2*s3 -> s2 => (0, -1, 0, 1)
s2*s3 -> s2*s1*s3*s2 => (1, 0, 0, -1)
s2*s3 -> s2*s3*s2 => (0, 0, 1, -1)
s2*s3 -> s3 => (0, -1, 1, 0)
s1*s2*s3 -> s1*s2 => (-1, 0, 0, 1)
s1*s2*s3 -> s1*s2*s1*s3*s2 => (0, 1, 0, -1)
s1*s2*s3 -> s1*s2*s3*s2 => (0, 0, 1, -1)
s1*s2*s3 -> s3 => (-1, 0, 1, 0)
s1*s2*s3 -> s2*s3 => (-1, 1, 0, 0)
Subgraph:
GKM graph with 6 nodes, valency 3 and axial function:
s2 -> id => (0, -1, 1, 0)
s1*s2 -> id => (-1, 0, 1, 0)
s1*s2 -> s2 => (-1, 1, 0, 0)
s3 -> id => (0, 0, -1, 1)
s2*s3 -> s2 => (0, -1, 0, 1)
s2*s3 -> s3 => (0, -1, 1, 0)
s1*s2*s3 -> s1*s2 => (-1, 0, 0, 1)
s1*s2*s3 -> s3 => (-1, 0, 1, 0)
s1*s2*s3 -> s2*s3 => (-1, 1, 0, 0)

julia> schubert_class(schubert, BO, "s1*s2")
(t3 - t4)*e[1] + (t2 - t4)*e[2] + (t1 - t4)*e[3]

julia> schubert_class(schubert, BO, "s1*s2*s3")
e[1] + e[2] + e[3] + e[4] + e[5] + e[6]

```
"""
function schubert_class(schubert::AbstractGKM_subgraph, BO::BruhatOrder, label::String)
  @req label in BO.labels "Label $label not found in Bruhat order."
  i = indexin([label], BO.labels)[1]
  return schubert_class(schubert, BO, i)
end


@doc raw"""
    schubert_class(schubert::AbstractGKM_graph, BO::BruhatOrder, v::Int64)

Like above, but the Weyl group element is given by the index `v` corresponding to its index as vertex of the GKM graph on the flag variety containing the Schubert variety.
"""
function schubert_class(schubert::AbstractGKM_subgraph, BO::BruhatOrder, v::Int64)
  descendants = Int64[]
  _add_all_descendants!(BO, v, descendants)
  subgraphVertices = [indexin([i], schubert.vDict)[1] for i in descendants]
  @req !any(a -> isnothing(a), subgraphVertices) "The demanded Schubert class is not contained in the given Schubert variety."
  return poincare_dual(gkm_subgraph_from_vertices(schubert.Domain, subgraphVertices))
end



"""
    schubert_class(flag::AbstractGKM_graph, BO::BruhatOrder, label::String)

Return the Poincare dual (as subvariety of the generalized partial flag variety `flag`) of the Schubert variety determined by `label`.

# Arguments
 - `flag::AbstractGKM_graph`: The generalized partial flag variety (as returned by `generalized_gkm_flag`).
- `BO::BruhatOrder`: The Bruhat order on the generalized partial flag variety (as returned by `get_bruhat_order_of_generalized_flag`).
- `label:String`: The label of the Weyl group element or coset defining the Schubert variety as subvariety of the generalized partial flag variety.

# Example
```jldoctest
julia> R = root_system(:G, 2);

julia> BO = get_bruhat_order_of_generalized_flag(R);

julia> flag = generalized_gkm_flag(R);

julia> schubert_class(flag, BO, "s2*s1*s2*s1*s2")
(-t2 + t3)*e[1] + (-t1 - t2 + 2*t3)*e[2] + (t1 - 2*t2 + t3)*e[3] + (-t1 + t3)*e[4] + (t1 - t2)*e[5] + (-t2 + t3)*e[7] + (-t1 - t2 + 2*t3)*e[8] + (t1 - 2*t2 + t3)*e[9] + (-t1 + t3)*e[10] + (t1 - t2)*e[11]

```
"""
function schubert_class(flag::AbstractGKM_graph, BO::BruhatOrder, label::String)
  @req label in BO.labels "Label $label not found in Bruhat order."
  i = indexin([label], BO.labels)[1]
  return schubert_class(flag, BO, i)
end

@doc raw"""
    schubert_class(flag::AbstractGKM_graph, BO::BruhatOrder, v::Int64)

Like above, but the Weyl group element is given by the index `v` corresponding to its index as vertex of the GKM graph on the Schubert variety.
"""
function schubert_class(flag::AbstractGKM_graph, BO::BruhatOrder, v::Int64)
  descendants = Int64[]
  _add_all_descendants!(BO, v, descendants)
  return poincare_dual(gkm_subgraph_from_vertices(flag, descendants))
end


@doc raw"""
    schubert_classes(schubert::AbstractGKM_subgraph, BO::BruhatOrder)

Return all Schubert classes on the given generalized Schubert variety.
The i-th row in the returned matrix is the Schubert class of the Weyl group element (or coset) corresponding
to the i-th vertex of the given Schubert variety.

# Arguments
- `schubert::AbstractGKM_subgraph`: The Schubert variety given as GKM subgraph of its generalized partial flag variety (as returned by `generalized_gkm_schubert`).
- `BO::BruhatOrder`: The Bruhat order for the generalized partial flag variety containing the Schubert variety (as returned by `get_bruhat_order_of_generalized_flag`).

# Examples
Schubert classes on the Schubert variety $\overline{X_{s_1s_2s_3}}\subset SL_4/P_1$:
```jldoctest schubert_classes
julia> R = root_system(:A, 2);

julia> S = simple_roots(R);

julia> schubert = generalized_gkm_schubert(R, S[1:1], "s1*s2");

julia> BO = get_bruhat_order_of_generalized_flag(R, S[1:1]);

julia> schubert_classes(schubert, BO)
3×3 Matrix{QQMPolyRingElem}:
 t1*t2 - t1*t3 - t2*t3 + t3^2  0        0
 t1 - t3                       t1 - t2  0
 1                             1        1
```
"""
function Oscar.schubert_classes(schubert::AbstractGKM_subgraph, BO::BruhatOrder)
  nv = n_vertices(schubert.Domain.g)
  z = zero(schubert.Domain.equivariantCohomology.coeffRing)
  M = fill(z, (nv, nv))
  for v in 1:nv
    sc = schubert_class(schubert, BO, schubert.Domain.labels[v])
    for j in 1:nv
      M[v, j] = sc[j]
    end
  end
  return M
end

@doc raw"""
    schubert_classes(flag::AbstractGKM_graph, BO::BruhatOrder)

Return all Schubert classes on the given generalized partial flag variety.
The i-th row in the returned matrix is the Schubert class of the Weyl group element (or coset) corresponding
to the i-th vertex of the given partial flag variety.

# Arguments
 - `flag::AbstractGKM_graph`: The generalized partial flag variety (as returned by `generalized_gkm_flag`).
 - `BO::BruhatOrder`: The Bruhat order for the generalized partial flag variety containing the Schubert variety (as returned by `get_bruhat_order_of_generalized_flag`).

# Example
```jldoctest schubert_classes_flag
julia> R = root_system(:A, 2);

julia> BO = get_bruhat_order_of_generalized_flag(R);

julia> flag = generalized_gkm_flag(R);

julia> schubert_classes(flag, BO)
6×6 Matrix{QQMPolyRingElem}:
 t1^2*t2 - t1^2*t3 - t1*t2^2 + t1*t3^2 + t2^2*t3 - t2*t3^2  0                             0        0  0                             0
 t1*t2 - t1*t3 - t2*t3 + t3^2                               t1*t2 - t1*t3 - t2*t3 + t3^2  0        0  0                             0
 t1 - t3                                                    t1 - t3                       t1 - t2  0  t1 - t2                       0
 1                                                          1                             1        1  1                             1
 t1^2 - t1*t2 - t1*t3 + t2*t3                               0                             0        0  t1^2 - t1*t2 - t1*t3 + t2*t3  0
 t1 - t3                                                    t2 - t3                       0        0  t1 - t3                       t2 - t3
```
"""
function Oscar.schubert_classes(flag::AbstractGKM_graph, BO::BruhatOrder)
  nv = n_vertices(flag.g)
  z = zero(flag.equivariantCohomology.coeffRing)
  M = fill(z, (nv, nv))
  for v in 1:nv
    sc = schubert_class(flag, BO, flag.labels[v])
    for j in 1:nv
      M[v, j] = sc[j]
    end
  end
  return M 
end

function _add_all_descendants!(BO::BruhatOrder, v::Int, descendants::Vector{Int64})
  (v in descendants) && return
  ord = BO.order
  for k in keys(ord)
    if k[1] == v
      push!(descendants, v)
      #println("Added $v")
      for d in ord[k]
        _add_all_descendants!(BO, d[1], descendants)
      end
    end
  end
end

function isrationally_smooth(base::AbstractGKM_graph, pt::String, rev::Bool=true)
    
  starting_elem = findfirst(v -> base.labels[v] == pt, 1:n_vertices(base.g))
  suborder = _create_order(base, starting_elem, rev, true)

  expected_valency = length(suborder[starting_elem, count('s', base.labels[starting_elem])])

  (expected_valency == count('s', pt)) || return false
#   println(suborder)

  for lab in keys(suborder)

    # number inbound
    inbound = count(t -> Base.in(lab, suborder[t]), keys(suborder))

    #number outbound
    outbound = length(suborder[lab])

    (inbound + outbound) != expected_valency && return false
    # println("$lab, $inbound $outbound")
  end

  return true
end
