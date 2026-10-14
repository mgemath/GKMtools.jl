struct _ConvexBundleFunctor end
struct _ConcaveBundleFunctor end

@inline _decorated_graph(dt::GW_decorated_tree) = dt.tree
@inline _decorated_graph(dt::GW_decorated_graph) = dt.g
@inline _edge_degrees(dt::GW_decorated_tree, e::Edge) = (edgeMult(e, dt),)
@inline _edge_degrees(dt::GW_decorated_graph, e::Edge) = edgeMult(e, dt)
@inline _decorated_valency(dt::GW_decorated_tree, v::Int) = degree(dt.tree, v)
@inline _decorated_valency(dt::GW_decorated_graph, v::Int) = valency(v, dt)

@inline function _bundle_top_weight(E::AbstractGKMVectorBundle, v::Int, t)
  result = one(t[1])
  @inbounds for i in 1:rank(E)
    result *= _bundle_weight_class(E, v, i, t)
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a >= 0 "virtual_zero_section is only implemented for convex vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in 0:(Int(a) * d)
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a < 0 "derivated_functor is only implemented for concave vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in (Int(a) * d + 1):-1
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result = result // _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result *= _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _gw_bundle_functor(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle, mode)
  @req baseof(E) == dt.gkm "The vector bundle and decorated graph must have the same base"
  if dt isa GW_decorated_graph
    @req all(iszero, dt.genus) "Only implemented for genus zero"
  end
  result = _bundle_edge_contribution(one(dt.context.t[1]), dt, E, mode)
  return _bundle_vertex_contribution(result, dt, E, mode)
end

"""
    virtual_zero_section(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the top Chern class of the genus-zero pushforward of `E`.
The bundle must be convex along every GKM edge.

# Examples
```jldoctest P7
julia> P7 = projective_space(NormalToricVariety, 7)
Normal toric variety
```

```jldoctest P7
julia> L = toric_line_bundle(P7, picard_group(P7)([2]))
Toric line bundle on a normal toric variety
```

```jldoctest P7
julia> l = gkm_vector_bundle_of_toric(L);

```

```jldoctest P7
julia> V = l + l + l + l;
```

```jldoctest P7
julia> X = baseof(V);
```

```jldoctest P7
julia> beta = curve_class(X, "1", "2");

julia> P = virtual_zero_section(V);

julia> gromov_witten(X, beta, 0, P; show_bar = false, fast_mode = true)
512

julia> gromov_witten(X, 2*beta, 0, P; show_bar = false, fast_mode = true)
9792
```

```jldoctest P3
julia> P3 = projective_space(NormalToricVariety, 3);

julia> O = trivial_line_bundle(P3);

julia> O3 = toric_line_bundle(P3, picard_group(P3)([3]));

julia> Fano = projectivization(O, O3); # Oscar's toric variety X

julia> L = anticanonical_bundle(Fano); # Oscar's anticanonical bundle

julia> V = gkm_vector_bundle_of_toric(L);

```

```jldoctest P3
julia> X = baseof(V);

julia> print_curve_classes(X)
2 -> 1: (0, 1), Chern number: 2
3 -> 1: (1, -3), Chern number: 1
4 -> 2: (1, 0), Chern number: 7
4 -> 3: (0, 1), Chern number: 2
5 -> 1: (1, -3), Chern number: 1
5 -> 3: (1, -3), Chern number: 1
6 -> 2: (1, 0), Chern number: 7
6 -> 4: (1, 0), Chern number: 7
6 -> 5: (0, 1), Chern number: 2
7 -> 1: (1, -3), Chern number: 1
7 -> 3: (1, -3), Chern number: 1
7 -> 5: (1, -3), Chern number: 1
8 -> 2: (1, 0), Chern number: 7
8 -> 4: (1, 0), Chern number: 7
8 -> 6: (1, 0), Chern number: 7
8 -> 7: (0, 1), Chern number: 2

julia> beta_1 = curve_class(X, "1", "2"); # first generator

julia> beta_2 = curve_class(X, "1", "3"); # second generator

julia> P = virtual_zero_section(V);

julia> gromov_witten(X, beta_1, 0, P; show_bar = false, fast_mode = true)
28

julia> gromov_witten(X, beta_2, 0, P; show_bar = false, fast_mode = true)
3
```

```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 1, 5);

julia> V = dual(S)^5;

julia> P4 = baseof(V);

julia> beta = curve_class(P4, "1", "2"); # line class

julia> P = virtual_zero_section(V);

julia> gromov_witten(P4, beta, 0, P; show_bar = false, fast_mode = true) # lines in the quintic in P4
2875

julia> gromov_witten(P4, 2*beta, 0, P; show_bar = false, fast_mode = true)
4876875//8

julia> gromov_witten(P4, 3*beta, 0, P; show_bar = false, fast_mode = true)
8564575000//27
```

```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 2, 4);

julia> plucker = dual(det(S));

julia> V = plucker^4;

julia> G24 = baseof(V);

julia> beta = curve_class(G24, "12", "13");

julia> P = virtual_zero_section(V);

julia> gromov_witten(G24, beta, 0, P; show_bar = false, fast_mode = true)
1280

julia> gromov_witten(G24, 2*beta, 0, P; show_bar = false, fast_mode = true)
92448

julia> gromov_witten(G24, 3*beta, 0, P; show_bar = false, fast_mode = true)
422690816//27
```

```jldoctest
julia> S, Q = tautological_and_univ_bd(GKMGraph, 2, 5);

julia> plucker = dual(det(S));

julia> V1 = plucker + (plucker^2) + (plucker^2);

julia> V2 = plucker + plucker + (plucker^3);

julia> V3 = dual(S) * plucker + (plucker^2);

julia> V4 = wedge_product(Q, 2) * plucker;

julia> P = virtual_zero_section.([V1, V2, V3, V4]);

julia> G25 = baseof(V1);

julia> beta = curve_class(G25, "12", "13");

julia> gromov_witten(G25, beta, 0, P; show_bar = false, fast_mode = true)
4-element Vector{QQFieldElem}:
 400
 540
 336
 325

julia> gromov_witten(G25, 2*beta, 0, P; show_bar = false, fast_mode = true)
4-element Vector{QQFieldElem}:
  5590
 25245//2
  3678
 25925//8
```
"""
function virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConvexBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConvexBundleFunctor()))
end

"""
    derivated_functor(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the top Chern class of the genus-zero first derived pushforward of `E`.
The bundle must be concave along every GKM edge.

# Examples
```jldoctest
julia> S, _ = tautological_and_univ_bd(GKMGraph, 1, 2); # S is the Serre's twisting bundle on P1

julia> V = S + S;

julia> P1 = baseof(V);

julia> beta = curve_class(P1, "1", "2"); # line class

julia> P = derivated_functor(V);

julia> gromov_witten(P1, beta, 0, P; show_bar = false, fast_mode = true)
1

julia> gromov_witten(P1, 2*beta, 0, P; show_bar = false, fast_mode = true)
1//8

julia> gromov_witten(P1, 3*beta, 0, P; show_bar = false, fast_mode = true)
1//27
```
"""
function derivated_functor(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConcaveBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConcaveBundleFunctor()))
end

"""
    reduced_virtual_zero_section(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the virtual zero section with the top Chern class at the last marked point divided out.

# Examples
```jldoctest
julia> P3 = projective_space(NormalToricVariety, 3);

julia> L = toric_line_bundle(P3, picard_group(P3)([2]));

julia> V = gkm_vector_bundle_of_toric(L);

julia> reduced_virtual_zero_section(V) isa GKMtools.EquivariantClass
true
```
"""
function reduced_virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_reduced_virtual_zero_section(dt, $E))
  return EquivariantClass(rule, dt -> _reduced_virtual_zero_section(dt, E))
end

function _reduced_virtual_zero_section(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle)
  @req !isempty(dt.marks) "Need at least one marked point to reduce the virtual zero section"
  result = _gw_bundle_functor(dt, E, _ConvexBundleFunctor())
  marked_vertex = imageOf(last(dt.marks), dt)
  return result // _bundle_top_weight(E, marked_vertex, dt.context.t)
end
