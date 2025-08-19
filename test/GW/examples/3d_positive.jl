# In this file, we compute some equivariant Seidel elements for the seven positive Hamiltonian 3D GKM graphs
# of complexity one that are not projections of toric examples (see [Charton--Kessler, Appendix A]).

# i is the position of the example in the list

E = Vector{Vector{QQBarFieldElem}}()

for i in 2:7

  println("3D positive non-toric no. $i:")

  G = gkm_3d_positive_non_toric(i)

  for j in 1:2
    println(Seidel_element(G, G.M[j], 2))
  end

  #print_curve_classes(G)
  #S = QH_structure_constants(G);

  #@req QH_is_polynomial(G) "QH for 3d non toric i=$i is not polynomial!"
  #@req QH_is_associative(G) "QH for 3d non toric i=$i is not associative"
  #@req QH_is_commutative(G) "QH for 3d non toric i=$i is not commutative"
  #@req QH_is_homogeneous(G) "QH for 3d non toric i=$i is not homogeneous"

  #push!(E, conjecture_O_eigenvalues(G))
  #println("Conjecture O eigenvalies for i=$i:")
  #for l in E[i]
  #  println("$l of absolute value $(abs(l))")
  #end

end