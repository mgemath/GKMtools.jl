export gromov_witten_pos_gen

function gromov_witten_pos_gen(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::EquivariantClass; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)
  return gromov_witten_pos_gen(G, beta, n_marks, max_genus, [P_input]; show_bar=show_bar, check_degrees=check_degrees, fast_mode)[1]
end

function gromov_witten_pos_gen(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::Array{EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)

  inputLength = length(P_input)

  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"
  @req max_genus >= 0 "Genus must be non-negative"
  @req n_marks >= 0 "Number of marks must be non-negative"

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

  max_edges::Int64 = _max_n_edges(H2, beta)

  # SHOW BAR WILL BE IMPLEMENTED LATER
#   if show_bar #set up progress data
#     number_trees = A000055(max_n_vert)
#     # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
#     # G.g, also counting ways to distribute the marked points among the vertices.
#     # This does not cound the edge multiplicities.
#     threshold = n_vertices(G.g) * sum(vert -> number_trees[vert] * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
#     progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
#     current_graph = 0
#  end

  # n_marks = length(classes)
  # iterate undecorated trees:

  for (genus, n_vert) in Iterators.product(0:max_genus, 2:(max_edges + 1)) # we fix the genus and the number of vertices
    
    (genus > 0) && n_vert < ceil(Int64, (3 + sqrt(1 + 8*genus)) / 2) && continue # skip impossible cases
    n_vert + genus - 1 > max_edges && continue # respect max_edges

    gen_array = Vector{Int64}[] # possible genus distributions on the vertices

    if max_genus - genus > 0
      for p in Combinatorics.partitions(max_genus - genus) 
          extended_p = vcat(p, zeros(Int64, n_vert - length(p)))
          for perm_extended_p in Combinatorics.multiset_permutations(extended_p, n_vert)
            push!(gen_array, perm_extended_p)
          end
      end
    else
      gen_array = [zeros(Int64, n_vert)]
    end

    for g6 in my_geng(genus, n_vert) # generation of graphs
        
      M = graph6_to_adjacency_matrix(g6)
      top_aut = compute_aut(M, n_vert) # automorphisms of the graph
      top_graph = graph_from_adjacency_matrix(Undirected, M) # graph in Oscar of the current iteration

      for col in Iterators.product([1:n_vertices(G.g) for _ in 1:nv(top_graph)]...) # iterate maps from graph to G.g:

        all(e -> col[src(e)] in nc[col[dst(e)]], edges(top_graph)) || continue # skip non-valid colorings

        Multi = [[1 for _ in 1:length(edges(top_graph))]] #_multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(top_graph)], beta) we do not need to iterate over edge multiplicities now, we do it later

        for gen_dist in gen_array # iterate genus distributions on the vertices

          for m_inv in Combinatorics.with_replacement_combinations(1:nv(top_graph), n_marks)  # iterate location of marks on the graph

            for edgeMult_array in Multi # iterate edge multiplicities

              PROD = prod(edgeMult_array)
              euler = zero(R.coeffRing)
              Euler = QQ(0)

              edgeMult = Dict{Edge, Int}(edges(top_graph) .=> edgeMult_array)
          
              for m in Combinatorics.multiset_permutations(m_inv, n_marks)

                ##### TEST
                println("Graph: $g6, aut:$(top_aut) Genus: $genus, n_vert: $n_vert, Coloring: $(collect(col)), Gen_dist: $gen_dist, Edge_mult: $edgeMult_array, Marks: $m")
                continue
                ##### END TEST

              #   dt = decoratedGraph(G, top_graph, collect(col), edgeMult, m)
            
              #   Class = [Base.invokelatest(P[k], dt) for k in keys(P)]

              #   all(c -> is_zero(c), Class) && continue

              #   if is_zero(euler) #euler == zero(R.coeffRing)
                  
              #     euler = Euler_inv(dt; check_degree=check_degrees)//(PROD * top_aut)
              #     for e in edges(top_graph)
              #       triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
              #       if !haskey(h_dict, triple)
              #         h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R; check=false, check_degrees=check_degrees)
              #       end
              #       euler *= h_dict[triple]
              #     end

              #     if fast_mode
              #       Euler = evaluate(euler, randoms)
              #     end

              #   end
              # end

              # if fast_mode
              #   foreach(i-> res[i] += evaluate(Class[i], randoms)*Euler, keys(Class)) 
              # else
              #   res += Class.*euler
              # end
              end
            end
          end
        end
      end
    end
  end
  return res
end
