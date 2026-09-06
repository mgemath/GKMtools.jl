@doc raw"""
    gromov_witten_nomarks(G, beta, classes; show_bar=true,
                          check_degrees=false, fast_mode=false,
                          threaded=false, g=0)
    gromov_witten_nomarks(G, beta, classes, n_marks, P_input;
                          show_bar=true, check_degrees=false, fast_mode=false,
                          threaded=false, g=0)

Integrate products of divisor-type insertions over the moduli space of
unmarked stable maps. Every entry of `classes` is a GKM class on `G` and is
integrated over the image curve of a stable map.

Passing a vector of vectors computes several products in one traversal of the
localization graphs. Positive genus is currently not supported by this
specialized formula.

The second form additionally multiplies the curve-integrated classes by
`P_input`, an `EquivariantClass` on a moduli space with `n_marks` marked
points.

Set `threaded=true` to distribute undecorated localization trees among the
available Julia threads. Each worker uses private result elements and caches.
The progress bar is disabled during threaded evaluation.

# Examples
julia> G24 = grassmannian(GKMGraph, 2, 4);

julia> beta = curve_class(G24, Edge(1, 2));

julia> gromov_witten_nomarks(G24, beta, [point_class(G24, 1), first_chern_class(G24), poincare_dual(subgraph_from_vertices(G24, [1, 2]))]; show_bar=false)
4

julia> P2 = projective_space(GKMGraph, 2);

julia> beta = curve_class(P2, Edge(1, 2));

julia> gromov_witten_nomarks(P2, beta, [point_class(P2, 1), point_class(P2, 1), point_class(P2, 1)]; show_bar=false)
2*t1 - t2 - t3

julia> gromov_witten_nomarks(P2, beta, [point_class(P2, 1)^2, point_class(P2, 1)]; show_bar=false)
t1^2 - t1*t2 - t1*t3 + t2*t3

julia> gromov_witten_nomarks(P2, beta, [point_class(P2, 1), point_class(P2, 1), point_class(P2, 1)]; show_bar=false)
2*t1 - t2 - t3

julia> gromov_witten_nomarks(P2, beta, [point_class(P2, 1), point_class(P2, 1), point_class(P2, 2)]; show_bar=false)
t1 - t3

julia> gromov_witten_nomarks(P2, beta, [point_class(P2, 1), point_class(P2, 1), point_class(P2, 3)]; show_bar=false)
t1 - t2
"""
function gromov_witten_nomarks(
  G::AbstractGKMGraph,
  beta::CurveClass,
  classes::AbstractVector{<:GKMClass};
  show_bar::Bool=true,
  check_degrees::Bool=false,
  fast_mode::Bool=false,
  threaded::Bool=false,
  g::Int64=0,
)
  return only(gromov_witten_nomarks(
    G, beta, [classes]; show_bar, check_degrees, fast_mode, threaded, g,
  ))
end

function _nomarks_result_degree(G, beta, classes)
  degree_cache = Tuple{GKMClass,Int}[]
  insertion_degree = sum(classes; init=0) do c
    cached_index = findfirst(entry -> first(entry) == c, degree_cache)
    if isnothing(cached_index)
      index = findfirst(!iszero, c.restrictions)
      degree = if isnothing(index)
        0
      else
        value = c.restrictions[index]
        _first_term_total_degree(numerator(value)) -
          _first_term_total_degree(denominator(value)) - 1
      end
      push!(degree_cache, (c, degree))
      degree
    else
      last(degree_cache[cached_index])
    end
  end
  virtual_dimension = valency(G) - 3 + Int(chern_number(G, beta))
  return insertion_degree - virtual_dimension
end

function _first_term_total_degree(f)
  iszero(f) && return 0
  ring = parent(f)
  ring in (ZZ, QQ) && return 0
  return sum(j -> exponent(f, 1, j), 1:nvars(ring); init=0)
end

function gromov_witten_nomarks(
  G::AbstractGKMGraph,
  beta::CurveClass,
  classes::AbstractVector{<:GKMClass},
  n_marks::Int64,
  P_input::EquivariantClass;
  show_bar::Bool=true,
  check_degrees::Bool=false,
  fast_mode::Bool=false,
  threaded::Bool=false,
  g::Int64=0,
)
  return only(gromov_witten_nomarks(
    G, beta, [classes], n_marks, P_input;
    show_bar, check_degrees, fast_mode, threaded, g,
  ))
end

function gromov_witten_nomarks(
  G::AbstractGKMGraph,
  beta::CurveClass,
  class_products::AbstractVector{<:AbstractVector{<:GKMClass}};
  show_bar::Bool=true,
  check_degrees::Bool=false,
  fast_mode::Bool=false,
  threaded::Bool=false,
  g::Int64=0,
)
  known_zero = fast_mode ? [
    _nomarks_result_degree(G, beta, classes) != 0 for classes in class_products
  ] : nothing
  return gromov_witten_nomarks(
    G, beta, class_products, 0, class_one();
    show_bar, check_degrees, fast_mode, threaded, g, known_zero,
  )
end

function gromov_witten_nomarks(
  G::AbstractGKMGraph,
  beta::CurveClass,
  class_products::AbstractVector{<:AbstractVector{<:GKMClass}},
  n_marks::Int64,
  P_input::Union{EquivariantClass,AbstractVector{<:EquivariantClass}};
  show_bar::Bool=true,
  check_degrees::Bool=false,
  fast_mode::Bool=false,
  threaded::Bool=false,
  g::Int64=0,
  known_zero::Union{Nothing,Vector{Bool}}=nothing,
)
  @req g >= 0 "Genus g must be non-negative."
  @req g == 0 "Positive genus is not yet supported for gromov_witten_nomarks."
  @req n_marks >= 0 "The number of marked points must be non-negative."
  @req !isempty(class_products) "gromov_witten_nomarks needs at least one input."
  if P_input isa AbstractVector
    @req length(P_input) == length(class_products) "Need one marked insertion per class product."
  end

  # Curve integrals depend only on the unoriented target edge. Cache them here
  # instead of recomputing them for every decorated tree.
  target_edges = collect(edges(graph(G)))
  edge_integral_cache = Tuple{GKMClass,Any,Any}[]
  edge_integral_cache_lock = ReentrantLock()
  insertions = map(enumerate(class_products)) do (input_index, product_classes)
    for c in product_classes
      @req graph(c) === G "All insertion classes must belong to G."
    end
    edge_integrals = fast_mode ? nothing : [
      _cached_nomarks_edge_integrals!(
        edge_integral_cache, edge_integral_cache_lock, c, G, target_edges, nothing,
      )
      for c in product_classes
    ]

    evaluate_insertion = function (dt)
      integrals_for_context = isnothing(edge_integrals) ? [
        _cached_nomarks_edge_integrals!(
          edge_integral_cache, edge_integral_cache_lock, c, G, target_edges, dt.context,
        ) for c in product_classes
      ] : edge_integrals
      prod(integrals_for_context; init=1) do integrals
        value = sum(edges(dt.tree); init=0) do e
          edgeMult(e, dt) * integrals[_unoriented_edge(imageOf(e, dt))]
        end
        value
      end
    end
    marked_insertion = P_input isa EquivariantClass ? P_input : P_input[input_index]
    EquivariantClass(:(_nomarks_insertion(dt)), evaluate_insertion) * marked_insertion
  end

  return _gromov_witten_gen_0(
    G, beta, n_marks, insertions, Val(fast_mode);
    show_bar, check_degrees, threaded, g, known_zero,
  )
end

@inline _unoriented_edge(e::Edge) =
  src(e) < dst(e) ? Edge(src(e), dst(e)) : Edge(dst(e), src(e))

function _cached_nomarks_edge_integrals!(cache, cache_lock, c, G, target_edges, context)
  return Base.lock(cache_lock) do
  # A vector cache deliberately uses == rather than object identity: callers
  # commonly construct the same point class several times in one expression.
  numeric_context = !isnothing(context) && eltype(context.t) === QQFieldElem
  cache_context = numeric_context ? context : nothing
  cached_index = findfirst(
    entry -> entry[1] == c && entry[2] === cache_context,
    cache,
  )
  isnothing(cached_index) || return cache[cached_index][3]

  integrals = numeric_context ?
    _specialized_nomarks_edge_integrals(c, G, target_edges, context) :
    _polynomial_nomarks_edge_integrals(c, G, target_edges)
  push!(cache, (c, cache_context, integrals))
  return integrals
  end
end

function _specialized_nomarks_edge_integrals(c, G, target_edges, context)
  values = _cached_restrictions!(context, G, c)
  result = Dict{Edge,QQFieldElem}()
  sizehint!(result, length(target_edges))
  for e in target_edges
    difference = values[src(e)] - values[dst(e)]
    result[_unoriented_edge(e)] = iszero(difference) ? zero(QQ) :
      difference / _weight_class(G, e, context.t)
  end
  return result
end

function _polynomial_nomarks_edge_integrals(c, G, target_edges)
  coefficient_ring = equivariant_coefficient_ring(G)
  t = collect(gens(coefficient_ring))
  values = c.restrictions

  # Polynomial GKM classes may be stored either directly or as fractions with
  # unit denominators. Exact division by the edge weight avoids costly
  # multivariate gcd normalization in the fraction field.
  polynomial_values = try
    map(values) do value
      parent(value) === coefficient_ring && return value
      denominator_value = denominator(value)
      is_unit(denominator_value) || throw(ArgumentError("localized restriction"))
      coefficient_ring(divexact(numerator(value), denominator_value))
    end
  catch
    nothing
  end

  if !isnothing(polynomial_values)
    result = Dict{Edge,eltype(polynomial_values)}()
    sizehint!(result, length(target_edges))
    for e in target_edges
      difference = polynomial_values[src(e)] - polynomial_values[dst(e)]
      result[_unoriented_edge(e)] = iszero(difference) ? zero(coefficient_ring) :
        divexact(difference, _weight_class(G, e, t))
    end
    return result
  end

  return Dict(_unoriented_edge(e) => integrate(c, e) for e in target_edges)
end
