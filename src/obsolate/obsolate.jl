function GKMproj_space(dim::Int; label::String = "x_")

  @warn "This feature is deprecated and may be removed in the future."
  
  g = complete_graph(dim+1)
  M = free_module(ZZ, dim+1)
  w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  for e in edges(g)
    w[e] = gens(M)[src(e)]-gens(M)[dst(e)]
  end
  labels = [label*"$i" for i in 0:dim]
  return gkm_graph(g, labels, M, w)
end
