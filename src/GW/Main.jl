@doc raw"""
    gromov_witten(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = false) -> GW invariants

Integrate the class `P_input` over the moduli space $\overline{\mathcal{M}_{0,n}}(X,\beta)$ of genus 0 stable maps to $X$ in class $\beta\in H_2(X;\mathbb{Z})$ with `n_marks`
marked points.
The result is an element of $\text{Frac}(H_T^*(\text{pt};\mathbb{Q}))$, i.e. a rational function in $\dim_\mathbb{C}(T)$ many variables.

!!! note
    If the underlying space is a (smooth projective or Hamiltonian) GKM space then the output should in fact live in $H_T^*(\text{pt};\mathbb{Q})$, so it should be 
    a polynomial in the $\dim_\mathbb{C}(T)$ many variables.

!!! warning
    The GKM graph `G` must have a connection, as this datum is required by the localization formula [LS17](@cite).

# Arguments
 - `G::AbstractGKM_graph`: The GKM graph of the target GKM space $X$.
 - `beta::CurveClass_type`: The (non-zero) curve class $\beta\in H_2(X;\mathbb{Z})$ in which the image of the stable map should lie.
    To produce `beta`, use functions like `curve_class` (see [Curve Classes](../GKM/CurveClasses.md)).
 - `P_input::EquivariantClass`: The equivariant cohomology class on $\overline{\mathcal{M}_{0,n}}(X,\beta)$ that is being integrated.
    Use the functions `ev`, `class_one`, and `Psi` to produce this. These classes also support arithmetic using `+`, `*`, et cetera.
 - `show_bar::Bool`: If `true`, a progress bar will be displayed showing the estimated time until completion. This should be used for big examples.
 - `fast_mode::Bool`: If the expected result of the computation is a number, this option will speed up the computation.

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
"""
function gromov_witten(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)
  return gromov_witten(G, beta, n_marks, [P_input]; show_bar=show_bar, check_degrees=check_degrees, fast_mode)[1]
end

function gromov_witten(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::Array{EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)

  inputLength = length(P_input)
  # inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology

  if fast_mode
    randoms = QQ.(rand(Int16, length(gens(R.coeffRing))))
    res = [zero(QQ) for _ in inputKeys] # zeros(QQFieldElem, inputSize)
  else
    res = [zero(R.coeffRingLocalized) for _ in inputKeys] # zeros(AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, inputSize)
  end
  
  
  if !is_effective(H2, beta)
    return res
  end

  P = [P_input[k].func for k in inputKeys]
  con = get_any_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  #########
  # Dict in order to store H
  h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
  ########

  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # G.g, also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = n_vertices(G.g) * sum(vert -> number_trees[vert] * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
 end

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(R.coeffRing)
          Euler = QQ(0)

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = [Base.invokelatest(P[k], dt) for k in keys(P)]

            all(c -> is_zero(c), Class) && continue

            if is_zero(euler) #euler == zero(R.coeffRing)
              euler = Euler_inv(dt; check_degree=check_degrees)//(PROD * aut)
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R; check=false, check_degrees=check_degrees)
                end
                euler *= h_dict[triple]
              end

              if fast_mode
                Euler = evaluate(euler, randoms)
              end

            end
            # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
            #@req _is_homogeneous(euler) "Euler not homogeneous"
            #@req _is_homogeneous(Class[1]) "Class not homogeneous"
            if fast_mode
              foreach(i-> res[i] += evaluate(Class[i], randoms)*Euler, keys(Class)) 
            else
              res += Class.*euler
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

@doc raw"""
    gromov_witten(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = false, check_degrees::Bool = false)

Same as before, but taking the total space of a GKM vector bundle as input.
"""
function gromov_witten(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = false, check_degrees::Bool = false)
  return gromov_witten(V, beta, n_marks, [P_input]; show_bar=show_bar, check_degrees=check_degrees)[1]
end

@doc raw"""
    gromov_witten(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass; show_bar::Bool = false, check_degrees::Bool = false)

Same as before, but taking the total space of a GKM vector bundle as input.
"""
function gromov_witten(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, P_input::Array{EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false)

  G = V.gkm

  inputLength = length(P_input)
  # inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input..."

  @req beta != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology
  res = [zero(R.coeffRingLocalized) for _ in inputKeys]
  if !is_effective(H2, beta)
    return res
  end

  P = [P_input[k].func for k in inputKeys]
  con = get_any_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"
  VBcon = get_connection(V)
  @req !isnothing(VBcon) "GKM vector bundle needs a connection!" #TODO: Could use any connection here.
  @req G.M == V.M "Gromov-Witten for vector bundles currently requires G.M == V.M"

  # Calculate weight classes in normal direction as elements of the coefficient ring
  _calculate_weight_classes(V)
  _calculate_connection_a(V)

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  #########
  # Dict in order to store H
  h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
  ########

  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # G.g, also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = n_vertices(G.g) * sum(vert -> number_trees[vert] * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
end

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(R.coeffRing)

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = [Base.invokelatest(P[k], dt) for k in inputKeys]
            #println("Tree($(dt.vDict)), Marks($(dt.marks)), Mult($(dt.edgeMult)):")
            #println(Class)
            all(c -> c==0, Class) && continue

            if euler == zero(R.coeffRing)
              euler = Euler_inv(dt; check_degree=check_degrees)//(PROD * aut)
              euler *= _Euler_inv_VB(dt, V)
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h_VB(V, Edge(col[src(e)], col[dst(e)]), triple[1], con, R; check=false, check_degrees=check_degrees)
                end
                euler *= h_dict[triple]
              end
            end

            res += Class.*euler
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