export store_colors, gromov_witten_load_color, gromov_witten_nomarks_load_color

function store_colors(G::AbstractGKM_graph, beta::CurveClass_type)

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  counter = 0
  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  name_file = "ls_color.txt"
  io = open(name_file, "a")

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    for col in CI   # colorings Iterator
      Serialization.serialize(io, (ls, col, parents, subgraph_ends))
      counter += 1
    end 
  end
  close(io)
  println("Stored $counter colorings to $name_file")
end


function gromov_witten_load_color(G::AbstractGKM_graph, beta::CurveClass_type, n_marks::Int64, P_input::Array{EquivariantClass}, line_num::Int64; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false)

  #########################
  ##  GENUS ZERO CASE:   ##
  #########################

  inputLength = length(P_input)
  # inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input."

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology

  if fast_mode
    res = [zero(QQ) for _ in inputKeys] # zeros(QQFieldElem, inputSize)

    # if fast_mode is activated, store edge weights and point euler classes locally.
    # These are passed to Euler_inv, _h, weight_class, and euler_class to optimize performance.
    edge_weight_dict = Dict{Edge, QQFieldElem}()
    point_weight_dict = vcat(Union{Nothing, QQFieldElem}[], repeat([nothing], n_vertices(G.g)))
    # t are the equivariant parameters.
    t = QQ.([(i+4)^2 for i in 1:length(gens(R.coeffRing))])

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

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########


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

  name_file = "ls_color.txt"
  pair = ([0], [0], [0], [0]) # dummy initialization
  io = open(name_file, "r")
  for i in 1:line_num
    if eof(io)
      error("File has fewer than $line_num lines.")
      return res
    end
    pair = deserialize(io)
    if i == line_num
      break
    end
  end
  close(io)

  (ls, col, parents, subgraph_ends) = pair
  # n_marks = length(classes)
  # iterate undecorated trees:
#   for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    # CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    # for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)
      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)
      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(t[1])

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = [Base.invokelatest(P[k], dt) for k in keys(P)]
            # TODO: can we pass t directly to each P[k]?

            all(c -> is_zero(c), Class) && continue

            # println("Class = $Class")

            if is_zero(euler) #euler == zero(R.coeffRing)
              euler = Euler_inv(dt, t, edge_weight_dict, point_weight_dict; check_degree=check_degrees)//(PROD * aut)
              #println("Euler: $euler")
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R, t, edge_weight_dict; check=false, check_degrees=check_degrees)
                end
                euler *= h_dict[triple]
              end

            end
            # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
            #@req _is_homogeneous(euler) "Euler not homogeneous"
            #@req _is_homogeneous(Class[1]) "Class not homogeneous"
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

        if show_bar #update the progress bar
          current_graph += tree_aut ÷ top_aut
          update!(progress_bar, current_graph,
              showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
        end
      end
#     end
#   end

  f = open("$line_num.txt", "w")
  println(f, res)
  close(f)
  return res
end

function gromov_witten_nomarks_load_color(G::AbstractGKM_graph, beta::CurveClass_type, evClasses::Vector{FreeModElem{QQMPolyRingElem}}, line_num::Int64; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false, g::Int64 = 0)
  res = gromov_witten_nomarks_load_color(G, beta, [evClasses], line_num; show_bar=show_bar, check_degrees=check_degrees, fast_mode=fast_mode, g=g)[1]
  namefile = res == 0 ? "$(line_num)_z.txt" : "$line_num.txt"
  f = open(namefile, "w")
  println(f, res)
  close(f)
  return res
end

function gromov_witten_nomarks_load_color(G::AbstractGKM_graph, beta::CurveClass_type, evClasses::Vector{Vector{FreeModElem{QQMPolyRingElem}}}, line_num::Int64; show_bar::Bool = true, check_degrees::Bool = false, fast_mode::Bool = false, g::Int64 = 0)

  @req g >= 0 "Genus g must be non-negative."
  # POSITIVE GENUS CASE: use functions in PosGen/Main_pos_gen.jl
  if g > 0
    @req g == 0 "Positive genus not yet supported for nomarks."
    #TODO: test in positive genus.
    return _gromov_witten_pos_gen(G, beta, n_marks, g, P_input; show_bar=show_bar, check_degrees = check_degrees, fast_mode=fast_mode)
  end

  #########################
  ##  GENUS ZERO CASE:   ##
  #########################

  @req !is_zero(beta) "Beta must be non-zero" # != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology
  
  inputLength = length(evClasses)
  inputKeys = keys(evClasses)
  if fast_mode
    res = [zero(QQ) for _ in inputKeys] # zeros(QQFieldElem, inputSize)

    # if fast_mode is activated, store edge weights and point euler classes locally.
    # These are passed to Euler_inv, _h, weight_class, and euler_class to optimize performance.
    edge_weight_dict = Dict{Edge, QQFieldElem}()
    point_weight_dict = vcat(Union{Nothing, QQFieldElem}[], repeat([nothing], n_vertices(G.g)))
    # t are the equivariant parameters.
    t = QQ.([(i+4)^2 for i in 1:length(gens(R.coeffRing))])
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

  con = get_any_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########


  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # G.g, also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = n_vertices(G.g) * sum(vert -> number_trees[vert] * ((length(nc[1]))^(vert - 1)) * binomial(vert+0-1, 0), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
 end
 
  name_file = "ls_color.txt"
  pair = ([0], [0], [0], [0]) # dummy initialization
  io = open(name_file, "r")
  for i in 1:line_num
    if eof(io)
      error("File has fewer than $line_num lines.")
      return res
    end
    pair = deserialize(io)
    if i == line_num
      break
    end
  end
  close(io)
  # n_marks = length(evClasses)
  # iterate undecorated trees:
  # for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
  
    (ls, col, parents, subgraph_ends) = pair
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    # CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    # for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), 0) # final zero is n_marks=0 for experimental version.
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(t[1])

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, 0))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = [prod(c -> _integrate(dt, c), k) for k in evClasses]
            # TODO: can we pass t directly to each P[k]?

            all(c -> is_zero(c), Class) && continue

            # println("Class = $Class")

            if is_zero(euler) #euler == zero(R.coeffRing)
              euler = Euler_inv(dt, t, edge_weight_dict, point_weight_dict; check_degree=check_degrees)//(PROD * aut)
              #println("Euler: $euler")
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R, t, edge_weight_dict; check=false, check_degrees=check_degrees)
                end
                euler *= h_dict[triple]
              end

            end
            # println("ls = $(ls), col = $(col), aut = $(aut), PRODW = $(PROD), m=$(m), E = $(euler)")
            #@req _is_homogeneous(euler) "Euler not homogeneous"
            #@req _is_homogeneous(Class[1]) "Class not homogeneous"
            if fast_mode
              # The isa(...) check below is necessary as sometimes Class[i] is an integer,
              # because evaluate(Int64, ...) is not defined.
              foreach(i-> res[i] += (isa(Class[i], Union{Number, QQFieldElem}) ? Class[i] : evaluate(Class[i], t))*euler, keys(Class)) 
              # foreach(i-> res[i] += (isa(Class[i], Union{Number, QQFieldElem}) ? Class[i] : evaluate(Class[i], t))*euler, keys(Class)) 
               # TODO: can we pass t directly to each P[k]? Then we don't need to evaluate here and get rid of this if-else block.
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
#     end
#   end

  return res
end