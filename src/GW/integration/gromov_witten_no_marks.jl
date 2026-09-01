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
  return gromov_witten_nomarks(
    G, beta, class_products, 0, class_one();
    show_bar, check_degrees, fast_mode, threaded, g,
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
  insertions = map(enumerate(class_products)) do (input_index, product_classes)
    edge_integrals = map(product_classes) do c
      @req graph(c) === G "All insertion classes must belong to G."
      Dict(_unoriented_edge(e) => integrate(c, e) for e in target_edges)
    end

    evaluate_insertion = function (dt)
      prod(edge_integrals; init=1) do integrals
        value = sum(edges(dt.tree); init=0) do e
          edgeMult(e, dt) * integrals[_unoriented_edge(imageOf(e, dt))]
        end
        _specialize_restriction(value, dt.context)
      end
    end
    marked_insertion = P_input isa EquivariantClass ? P_input : P_input[input_index]
    EquivariantClass(:(_nomarks_insertion(dt)), evaluate_insertion) * marked_insertion
  end

  return gromov_witten(
    G, beta, n_marks, insertions;
    show_bar, check_degrees, fast_mode, threaded, g,
  )
end

@inline _unoriented_edge(e::Edge) =
  src(e) < dst(e) ? Edge(src(e), dst(e)) : Edge(dst(e), src(e))
