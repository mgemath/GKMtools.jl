@doc raw"""
    gromov_witten(G::AbstractGKMGraph, beta::CurveClass, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = true, fast_mode::Bool = false, g::Int64 = 0) -> GW invariants

Integrate the class `P_input` over the moduli space $\overline{\mathcal{M}}_{g,n}(X,\beta)$ of genus `g` stable maps to $X$ in class $\beta\in H_2(X;\mathbb{Z})$ with `n_marks`
marked points.
The result is an element of $\text{Frac}(H_T^*(\text{pt};\mathbb{Q}))$, i.e. a rational function in $\dim_\mathbb{C}(T)$ many variables.

!!! note
    If the underlying space is a (smooth projective or Hamiltonian) GKM space then the output should in fact live in $H_T^*(\text{pt};\mathbb{Q})$, so it should be
    a polynomial in the $\dim_\mathbb{C}(T)$ many variables.

!!! warning
    The GKM graph `G` must have a connection, as this datum is required by the localization formula [LS17](@cite).

# Arguments
 - `G::AbstractGKMGraph`: The GKM graph of the target GKM space $X$.
 - `beta::CurveClass`: The (non-zero) curve class $\beta\in H_2(X;\mathbb{Z})$ in which the image of the stable map should lie.
    To produce `beta`, use functions like `curve_class` (see [Curve Classes](../GKM/CurveClasses.md)).
 - `P_input::EquivariantClass`: The equivariant cohomology class on $\overline{\mathcal{M}}_{g,n}(X,\beta)$ that is being integrated.
    Use the functions `ev`, `class_one`, and `Psi` to produce this. These classes also support arithmetic using `+`, `*`, et cetera.
 - `show_bar::Bool`: If `true`, a progress bar will be displayed showing the estimated time until completion. This should be used for big examples.
 - `fast_mode::Bool`: If the expected result of the computation is a number, this option will speed up the computation.
 - `g::Int64`: Genus of the GW invariant to be computed. The default value is `0`.

!!! warning
    If the expected result of the computation is not a number and `fast_mode` is `true`, the result will be a meaningless number.

# Example
```jldoctest gromov_witten
julia> P2 = projective_space(GKM_graph, 2);

julia> beta = curve_class(P2, Edge(1, 2));

julia> gromov_witten(P2, beta, 1, ev(1, point_class(P2, 1)); show_bar=false)
0

julia> gromov_witten(P2, beta, 2, ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 2)); show_bar=false)
1

julia> gromov_witten(P2, beta, 2, ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 1)); show_bar=false, fast_mode=true)
1

julia> gromov_witten(P2, beta, 2, ev(1, point_class(P2, 1))^2 * ev(2, point_class(P2, 2)); show_bar=false)
t1^2 - t1*t2 - t1*t3 + t2*t3

julia> gromov_witten(P2, beta, 3, ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 1)) * ev(3, point_class(P2, 3)); show_bar=false)
t1 - t2
```

# Example in positive genus

In the following example, $F$ is the twisted flag manifold and $\beta=[C_e]$, where
$e$ is the unique edge of the GKM graph of $F$ with $\int_{C_e} c_1(T_F)=0$.
Thus, the resulting invariants agree with those of an equivariantly Calabi--Yau GKM linearization of
$\mathcal{O}_{\mathbb{P}^1}(1)\oplus \mathcal{O}_{\mathbb{P}^1}(-3)$.

```jldoctest gromov_witten
julia> F = gkm_3d_twisted_flag();

julia> beta = curve_class(F, "3", "4")
(-2, 1)

julia> for d in 1:3
         gw = gromov_witten(F, d*beta, 0, class_one(); g=1, show_bar=false)
         println("Genus 1, degree $d: $gw")
       end
Genus 1, degree 1: 1//12
Genus 1, degree 2: -1//24
Genus 1, degree 3: -29//36
```
"""
function gromov_witten(G::AbstractGKMGraph, beta::CurveClass, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false, g::Int64 = 0)
  return gromov_witten(G, beta, n_marks, [P_input]; show_bar=show_bar, check_degrees=check_degrees, fast_mode, g=g)[1]
end

function gromov_witten(G::AbstractGKMGraph, beta::CurveClass, n_marks::Int64, P_input::AbstractVector{<:EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false, g::Int64 = 0)
  @req g >= 0 "Genus g must be non-negative."
  g > 0 && return _gromov_witten_pos_gen(G, beta, n_marks, g, P_input; show_bar, check_degrees, fast_mode)
  return _gromov_witten_gen_0(G, beta, n_marks, P_input, Val(fast_mode); show_bar, check_degrees, g)
end

function _gromov_witten_gen_0(G::AbstractGKMGraph, beta::CurveClass, n_marks::Int64, P_input::AbstractVector{<:EquivariantClass}, ::Val{fast_mode}; show_bar::Bool, check_degrees::Bool, g::Int64) where fast_mode

  @req g >= 0 "Genus g must be non-negative."
  # POSITIVE GENUS CASE: use functions in PosGen/Main_pos_gen.jl
  # if g > 0
  #   return _gromov_witten_pos_gen(G, beta, n_marks, g, P_input; show_bar=show_bar, check_degrees = check_degrees, fast_mode=fast_mode)
  # end

  #########################
  ##  GENUS ZERO CASE:   ##
  #########################

  inputLength = length(P_input)
  # inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = get_cohomology(G)

  if fast_mode
    res = [zero(QQ) for _ in inputKeys] # zeros(QQFieldElem, inputSize)

    # t are the equivariant parameters.
    t = QQ.(rand(Int16, length(gens(R.coefficient_ring))))

    #########
    # Dict in order to store H
    h_dict = Dict{Tuple{Int64, Int64, Int64}, QQFieldElem}() # Lambda_gamma_e_dict
    class_context = GWClassEvaluationContext(t, QQFieldElem)
    ########
  else
    res = [zero(R.localized_coefficient_ring) for _ in inputKeys] # zeros(AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, inputSize)

    # t are the equivariant parameters.
    t = gens(R.coefficient_ring)

    #########
    # Dict in order to store H
    h_dict = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
    class_context = GWClassEvaluationContext(t, typeof(zero(R.localized_coefficient_ring)))
    ########
  end


  if !is_effective(H2, beta)
    return res
  end

  P = map(ec -> ec.func, P_input)
  con = connection(G)
  multiplicity_cache = Dict{Tuple{Vararg{Edge}}, Set{Vector{Int}}}()
  # @req !isnothing(con) "GKM graph needs a connection!"

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(graph(G))
    nc[v] = sort(all_neighbors(graph(G), v))
  end
  #########


  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # graph(G), also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = n_vertices(graph(G)) * sum(vert -> number_trees[vert] * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
 end

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_edges = collect(edges(tree))
    tree_aut = count_iso(ls)

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to graph(G):
    for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      target_edges = [Edge(col[src(e)], col[dst(e)]) for e in tree_edges]
      multiplicity_key = Tuple(Edge(min(src(e), dst(e)), max(src(e), dst(e))) for e in target_edges)
      Multi = get!(multiplicity_cache, multiplicity_key) do
        _multiplicities(H2, target_edges, beta)
      end

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)

        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(t[1])
          euler_computed = false

          edgeMult = Dict{Edge, Int}(tree_edges .=> edgeMult_array)

          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))


            dt = decoratedTree(G, tree, col, edgeMult, m, class_context; check=false)

            Class = [P[k](dt) for k in eachindex(P)]

            all(c -> is_zero(c), Class) && continue

            # println("Class = $Class")

            if !euler_computed
              euler = Euler_inv(dt, t; check_degree=check_degrees)//(PROD * aut)
              #println("Euler: $euler")
              for e in tree_edges
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                h = get!(h_dict, triple) do
                  _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, G, t)
                end
                euler *= h
              end
              euler_computed = true

            end
            # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
            #@req _is_homogeneous(euler) "Euler not homogeneous"
            #@req _is_homogeneous(Class[1]) "Class not homogeneous"
            if fast_mode
              # The isa(...) check below is necessary as sometimes Class[i] is an integer,
              # because evaluate(Int64, ...) is not defined.
              for i in eachindex(Class)
                value = Class[i]
                res[i] += (value isa Union{Number, QQFieldElem} ? value : evaluate(value, t)) * euler
              end
              else
              for i in eachindex(Class)
                res[i] += Class[i] * euler
              end
            end

          end
        end

        if show_bar #update the progress bar
          current_graph += tree_aut ÷ top_aut
          update!(progress_bar, current_graph,
              showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
        end
      end
    end
  end
  return res
end

# For debugging purposes:
function _get_degree(f)
  if f == 0
    return nothing
  end
  f = f//1
  return _get_deg_poly(numerator(f)) - _get_deg_poly(denominator(f))
end

function _is_polynomial(f)
  if f == 0
    return true
  end
  f = f//1
  return _get_deg_poly(denominator(f)) == 0
end

function _get_deg_poly(f)
  for e in exponents(f)
    return sum(e)
  end
end

function _is_homogeneous(f)
  f = f//1
  return _is_homogeneous_poly(numerator(f)) && _is_homogeneous_poly(denominator(f))
end

function _is_homogeneous_poly(f)
  if f == 0
    return true
  end
  s::Union{Nothing, Int64} = nothing
  for e in exponents(f)
    if isnothing(s)
      s = sum(e)
    else
      if s != sum(e)
        return false
      end
    end
  end
  return true
end
