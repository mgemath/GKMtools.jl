export _gromov_witten_pos_gen

function _gromov_witten_pos_gen(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::EquivariantClass; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)
  return _gromov_witten_pos_gen(G, beta, n_marks, max_genus, [P_input]; show_bar=show_bar, check_degrees=check_degrees, fast_mode=fast_mode)[1]
end

function _gromov_witten_pos_gen(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::Array{EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)

  inputLength = length(P_input)

  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"
  @req max_genus >= 0 "Genus must be non-negative"
  @req n_marks >= 0 "Number of marks must be non-negative"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology

  if fast_mode
    res = [zero(QQ) for _ in inputKeys] # zeros(QQFieldElem, inputSize)

    # if fast_mode is activated, store edge weights and point euler classes locally.
    # These are passed to Euler_inv, _h, weight_class, and euler_class to optimize performance.
    edge_weight_dict = Dict{Edge, QQFieldElem}()
    point_weight_dict = vcat(Union{Nothing, QQFieldElem}[], repeat([nothing], n_vertices(G.g)))
    # TODO (D): change to Giosue's version (or use zero instead of nothing.)
    # t are the equivariant parameters.
    t = QQ.(rand(Int16, length(gens(R.coeffRing))))

    #########
    # Dict in order to store H
    h_dict = Dict{Tuple{Int64, Int64, Int64}, QQFieldElem}() # Lambda_gamma_e_dict
    ########
  else
    res = [zero(R.coeffRingLocalized) for _ in inputKeys] # zeros(AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, inputSize)

    # if we are not in fast_mode, edge weights and point euler classes are polynomials in the equivariant parameters,
    # which are already stored in R.
    edge_weight_dict = R.edgeWeightClasses
    point_weight_dict = R.pointEulerClasses
    # t are the equivariant parameters.
    t = gens(R.coeffRing)

    #########
    # Dict in order to store H
    h_dict = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
    ########
  end


  if !is_effective(H2, beta)
    return res
  end

  P = [P_input[k].func for k in inputKeys]
  con = get_any_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"

  ##### check psi classes in P_input
  sum_psi_constant = max([P_input[k].sum_psi for k in inputKeys]...)
  prod_psi_constant = max([P_input[k].prod_psi for k in inputKeys]...)
  @req sum_psi_constant == false "In positive genus you cannot sum two equivariant classes where psi appears just in one"
  @req prod_psi_constant == false "You cannot multiply two equivariant classes where psi appears in both"

  psi_exp = P_input[inputKeys[1]].psi_exp ;
  lambda_coef = P_input[inputKeys[1]].lambda_coef ;

  @req all(k -> P_input[k].psi_exp == psi_exp, inputKeys) "All equivariant classes must have the same psi exponents in positive genus."
  @req all(k -> P_input[k].lambda_coef == lambda_coef, inputKeys) "All equivariant classes must have the same lambda coefficients."
  
  @req length(lambda_coef) <= max_genus "We must have length(lambda_coef) <= max_genus"
  @req length(psi_exp) <= n_marks "We must have length(psi_exp) <= n_marks"
  #####

  ##### add zeros to psi_exp and lambda_coef to ensure correct length
  lambda_coef = vcat(lambda_coef, zeros(Int64, max_genus - length(lambda_coef)))
  psi_exp = vcat(psi_exp, zeros(Int64, n_marks - length(psi_exp)))
  have_psis = any(a -> !iszero(a), psi_exp)
  println("lambda_coef = $lambda_coef")
  println("psi_exp = $psi_exp")
  #####

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  max_edges::Int64 = _max_n_edges(H2, beta)
  # VPs = vertex_polynomials(max_genus, max_edges, valency(G))
  H = load_H(max_genus, max_edges + n_marks)

  ## Progress bar
  if show_bar
    threshold = compute_threshold_for_progress_bar(max_genus, max_edges, nc)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
  end

  for (top_genus, n_vert) in Iterators.product(0:max_genus, 2:(max_edges + 1)) # we fix the top_genus and the number of vertices

    n_vert + top_genus - 1 > max_edges && continue # respect max_edges
    (top_genus > 0) && n_vert < ceil(Int64, (3 + sqrt(1 + 8*top_genus)) / 2) && continue # skip impossible cases

    for g6 in my_geng(top_genus, n_vert) # generation of graphs

      M = graph6_to_adjacency_matrix(g6)
      top_graph = graph_from_adjacency_matrix(Undirected, M) # graph in Oscar of the current iteration
      top_graph_Graphs = Graphs.SimpleGraph(M) # graph in Graphs of the current iteration

      top_aut = compute_aut(top_graph_Graphs) # automorphisms of the graph

      for (col, col_aut) in collect_all_uniques_cols_W_C_N_flatmap(top_graph_Graphs, nc, top_aut) # iterate colorings modulo isomorphisms, return a pair (coloring, automorphisms of the coloring)

        Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(top_graph)], beta)

        for (gen_dist, gen_dist_aut) in Iterators.flatmap(multiedge_grow -> genus_distribution_mod_iso(top_graph_Graphs, col, col_aut, max_genus, top_genus + multiedge_grow), 0:(max_genus - top_genus)) # iterate genus distributions on the vertices

          for (multiedges, multiedges_aut) in multiedges_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, gen_dist_aut, max_genus, top_genus, Multi) # TODO: what is multiedges_aut?          
            
            PROD = prod(prod.(multiedges))
            aut = compute_internal_aut_multiedges(multiedges) * multiedges_aut
            edgeMult = Dict{Edge, Vector{Int64}}(edges(top_graph) .=> multiedges)

            for m_inv in Combinatorics.with_replacement_combinations(1:nv(top_graph), n_marks)  # iterate location of marks on the graph

              euler = zero(t[1])

              for m in Combinatorics.multiset_permutations(m_inv, length(m_inv))

                
                # println("Graph: $g6, aut:$multiedges_aut Genus: $top_genus, Coloring: $col, Gen_dist: $gen_dist, Edge_mult: $multiedges, PROD=$PROD")
                # println("Marks: $m")
                # println("Total aut: $aut")
                # continue

                dg = decoratedGraph(G, top_graph, col, edgeMult, m, gen_dist)

                Class = [Base.invokelatest(P[k], dg) for k in keys(P)]

                # all(c -> is_zero(c), Class) && continue

                if have_psis || is_zero(euler) # If we have psis, then relabeling the marked points gives different Euler_inv_pos_gen.

                  euler = Euler_inv_pos_gen(dg, t, edge_weight_dict, point_weight_dict, psi_exp, H)//(PROD * aut)
                  # println("Euler (w/o h) = $(factor(numerator(euler))) // $(factor(denominator(euler)))")
                  for e in edges(top_graph)
                    for em in dg.edgeMult[e]
                      triple = (em, min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                      if !haskey(h_dict, triple)
                        h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R, t, edge_weight_dict; check=false, check_degrees=check_degrees)
                      end
                      euler *= h_dict[triple]
                      # println("h = $(factor(numerator(h_dict[triple]))) // $(factor(denominator(h_dict[triple])))")
                    end
                  end
                end

                if fast_mode
                  # The isa(...) check below is necessary as sometimes Class[i] is an integer,
                  # because evaluate(Int64, ...) is not defined.
                  foreach(i-> res[i] += (isa(Class[i], Union{Number, QQFieldElem}) ? Class[i] : evaluate(Class[i], t))*euler, keys(Class))
                  # TODO: can we pass t directly to each P[k]? Then we don't need to evaluate here and get rid of this if-else block.
                else
                  res += Class.*euler
                end

              end
            end
          
          end

          if show_bar #update the progress bar
            current_graph += top_aut ÷ gen_dist_aut #col_aut
            update!(progress_bar, current_graph,
              showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
          end

        end

        
      end
    end
  end
  return res
end


# function _gromov_witten_pos_gen(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::EquivariantClass; show_bar::Bool = false, check_degrees::Bool = false)
#   return _gromov_witten_pos_gen(V, beta, n_marks, max_genus, [P_input]; show_bar=show_bar, check_degrees=check_degrees)[1]
# end
# 
# function _gromov_witten_pos_gen(V::GKM_vector_bundle, beta::CurveClass_type, n_marks::Int64, max_genus::Int64, P_input::Array{EquivariantClass}; show_bar::Bool = true, check_degrees::Bool = false)
# 
#   G = V.gkm
#   R = G.equivariantCohomology
# 
#   edge_weight_dict = R.edgeWeightClasses
#   point_weight_dict = R.pointEulerClasses
#   t = gens(R.coeffRing)
# 
#   inputLength = length(P_input)
#   # inputSize = size(P_input)
#   inputKeys = keys(P_input)
#   @req inputLength > 0 "gromov_witten needs at least one input for P_input..."
# 
#   @req beta != zero(parent(beta)) "Beta must be non-zero"
# 
#   H2 = GKM_second_homology(G)
#   R = G.equivariantCohomology
#   res = [zero(R.coeffRingLocalized) for _ in inputKeys]
#   if !is_effective(H2, beta)
#     return res
#   end
# 
#   P = [P_input[k].func for k in inputKeys]
#   con = get_any_connection(G)
#   @req !isnothing(con) "GKM graph needs a connection!"
#   VBcon = get_any_connection(V)
#   @req !isnothing(VBcon) "GKM vector bundle needs a connection!"
#   @req G.M == V.M "Gromov-Witten for vector bundles currently requires G.M == V.M"
# 
#   # Calculate weight classes in normal direction as elements of the coefficient ring
#   _calculate_weight_classes(V)
#   _calculate_connection_a(V)
# 
#   ########
#   # this part is needed for the generation of colorings
#   nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
#   for v in 1:n_vertices(G.g)
#     nc[v] = sort(all_neighbors(G.g, v))
#   end
#   #########
# 
#   #########
#   # Dict in order to store H
#   h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
#   ########
# 
#   max_edges::Int64 = _max_n_edges(H2, beta)
#   VPs = vertex_polynomials(max_genus, max_edges, valency(G) + rank(V))
# 
#   ## Progress bar
#   if show_bar
#     threshold = compute_threshold_for_progress_bar(max_genus, max_edges)
#     progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
#     current_graph = 0
#   end
# 
#   for (top_genus, n_vert) in Iterators.product(0:max_genus, 2:(max_edges + 1)) # we fix the top_genus and the number of vertices
# 
#     n_vert + top_genus - 1 > max_edges && continue # respect max_edges
#     (top_genus > 0) && n_vert < ceil(Int64, (3 + sqrt(1 + 8*top_genus)) / 2) && continue # skip impossible cases
# 
#     for g6 in my_geng(top_genus, n_vert) # generation of graphs
# 
#       M = graph6_to_adjacency_matrix(g6)
#       top_graph = graph_from_adjacency_matrix(Undirected, M) # graph in Oscar of the current iteration
#       top_graph_Graphs = Graphs.SimpleGraph(M) # graph in Graphs of the current iteration
# 
#       top_aut = compute_aut(top_graph_Graphs) # automorphisms of the graph
# 
#       for (gen_dist, gen_dist_aut) in Iterators.flatmap(multiedge_grow -> genus_distribution_mod_iso(top_graph_Graphs, top_aut, max_genus, top_genus + multiedge_grow), 0:(max_genus - top_genus)) # iterate genus distributions on the vertices
# 
#         for (col, col_aut) in colorings_modulo_iso(top_graph_Graphs, nc, gen_dist, gen_dist_aut) # iterate colorings modulo isomorphisms, return a pair (coloring, automorphisms of the coloring)
# 
#           for (multiedges, multiedges_aut) in multiedges_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, col_aut, max_genus, top_genus, beta, H2) # TODO: what is multiedges_aut?
# 
#             PROD = prod(prod.(multiedges))
#             aut = compute_internal_aut_multiedges(multiedges) * multiedges_aut
#             edgeMult = Dict{Edge, Vector{Int64}}(edges(top_graph) .=> multiedges)
# 
#             for m_inv in Combinatorics.with_replacement_combinations(1:nv(top_graph), n_marks)  # iterate location of marks on the graph
# 
#               euler = zero(t[1])
# 
#               for m in Combinatorics.multiset_permutations(m_inv, length(m_inv))
# 
# 
#                 # println("Graph: $g6, aut:$multiedges_aut Genus: $top_genus, Coloring: $col, Gen_dist: $gen_dist, Edge_mult: $multiedges, PROD=$PROD")
#                 # println("Marks: $m")
#                 # println("Total aut: $aut")
#                 # continue
# 
#                 dg = decoratedGraph(G, top_graph, col, edgeMult, m, gen_dist)
# 
#                 Class = [Base.invokelatest(P[k], dg) for k in keys(P)]
# 
#                 all(c -> is_zero(c), Class) && continue
# 
#                 if is_zero(euler)
# 
#                   euler = _Euler_inv_pos_gen_VB(dg, V, VPs)//(PROD * aut)
#                   # println("Euler (w/o h) = $(factor(numerator(euler))) // $(factor(denominator(euler)))")
#                   for e in edges(top_graph)
#                     for em in dg.edgeMult[e]
#                       triple = (em, min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
#                       if !haskey(h_dict, triple)
#                         h_dict[triple] = _h_VB(V, Edge(col[src(e)], col[dst(e)]), triple[1], con, R; check=false, check_degrees=check_degrees)
#                       end
#                       euler *= h_dict[triple]
#                       # println("h = $(factor(numerator(h_dict[triple]))) // $(factor(denominator(h_dict[triple])))")
#                     end
#                   end
#                 end
# 
#                 res += Class.*euler
#               end
#             end
#           end
#         end
# 
#         if show_bar #update the progress bar
#           current_graph += top_aut ÷ gen_dist_aut
#           update!(progress_bar, current_graph,
#               showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
#         end
#       end
#     end
#   end
#   return res
# 
# end