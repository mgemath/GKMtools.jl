export gromov_witten_part

function gromov_witten_part(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::Array{EquivariantClass}, col_min::Int64, col_max::Int64, scal::Vector{Int64}=Int64[]; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)

  inputLength = length(P_input)
  # inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"

  @req col_min <= col_max "col_min must be less or equal to col_max"
  @req col_min >= 1 "col_min must be at least 1"
  @req col_max <= n_vertices(G.g) "col_max must be at most the number of vertices of the GKM graph"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology

  if fast_mode
    randoms = isempty(scal) ? QQ.([(i+4)^2 for i in 1:length(gens(R.coeffRing))]) : QQ.(scal)
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

  println("Computing Gromov-Witten invariant between vertices $col_min and $col_max.")

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

            if col[1] < col_min || col[1] > col_max
                Multi = Vector{Vector{Int64}}()
            else
                Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)
            end

            # Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

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

function gromov_witten_part(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::EquivariantClass, col_min::Int64, col_max::Int64, scal::Vector{Int64}=Int64[]; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)
  return gromov_witten_part(G, beta, n_marks, [P_input], col_min, col_max, scal; show_bar=show_bar, check_degrees=check_degrees, fast_mode)[1]
end
