#TODO do we still need this now that we have _multiplicities() in curveClasses.jl?

function all_ones_multip(G::AbstractGKM_graph, tree::Graph)

    return [Dict{Edge, Int}(edges(tree) .=> [1 for _ in 1:n_edges(tree)])]
end

function multi(G::AbstractGKM_graph, tree::Graph, beta)
  if bettiNumbers(G)[2] == 1 && typeof(beta) == Int64
    return collect(Iterators.Flatten([Combinatorics.multiset_permutations(p, n_edges(tree)) for p in Combinatorics.partitions(beta, n_edges(tree))]))
  end

  return nothing
end