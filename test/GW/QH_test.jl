G = projective_space(GKM_graph, 3)
#G = blowupGKM(GKMsubgraph_from_vertices(G, [1, 2])).super
beta = curve_class(G, Edge(1, 2))
class1 = point_class(1, G)
class2 = point_class(2, G)
#class3 = pointClass(3, G)
t1, t2 = gens(G.equivariantCohomology.coeffRing)

a = quantumProduct(G, beta, class1, class2)
b = quantumProduct(G, beta, class1, class2; useStructureConstants=false)
#println(a)
#println(b)
println(a==b)
a = quantumProduct(G, beta, class1*class1, one(G.equivariantCohomology))
b = quantumProduct(G, beta, class1*class1, one(G.equivariantCohomology); useStructureConstants=false)
println(a==b)

S = GKMtools.QH_structure_constants(G)
Sfactored = Dict{GKMtools.CurveClass_type, Any}()
for beta in keys(S)
  println("Structure constants for beta=$beta:")
  Sfactored[beta] = [S[beta][k] == 0 ? 0 : factor(numerator(S[beta][k])) for k in keys(S[beta])]
  println(Sfactored[beta])
end

# P1 example in nice base:

P1 = GKMproj_space(1)
beta1 = curve_class(P1, Edge(1, 2))
u0, u1 = gens(P1.equivariantCohomology.coeffRing)
S1 = GKMtools.QH_structure_constants(P1)

base1 = [1 1; u0-u1 0 ]
SB1 = GKMtools.QH_structure_constants_in_basis(P1, base1)

# P2 example in nice base:

P2 = GKMproj_space(2)
beta = curve_class(P2, Edge(1, 2))
u0, u1, u2 = gens(P2.equivariantCohomology.coeffRing)
S2 = GKMtools.QH_structure_constants(P2)

base = [1 1 1 ; u0-u2 u1-u2 0 ; (u0-u1)*(u0-u2) 0 0 ]
SB2 = GKMtools.QH_structure_constants_in_basis(P2, base)

H = [u0, u1, u2]
quantumProduct(P2, beta, H.-u1, (H.-u2).*(H.-u0)) # evaluates to one, as required!
#New notation for this:
H = QH_class(P2, H)
println((H-u0)*(H-u1)*(H-u2))

# Check conjecture O:
println(conjecture_O_eigenvalues(P2))