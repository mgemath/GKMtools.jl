# In this file we study the smooth Schubert variety in type A3 in G/P given by w=s1*s3*s2 and P given by a1.

# To see this Schubert variety as subvariety of Fl([2,1,1]), use:

# julia> F211 = flag_variety(GKM_graph, [2,1,1]);

# julia> S = gkm_subgraph_from_vertices(F, ["123", "132", "124", "231", "142", "241"]).self
# GKM graph with 6 nodes, valency 3 and axial function:
# 124 -> 123 => (0, 0, -1, 1)
# 132 -> 123 => (0, -1, 1, 0)
# 142 -> 124 => (0, -1, 0, 1)
# 142 -> 132 => (0, 0, -1, 1)
# 231 -> 123 => (-1, 0, 1, 0)
# 231 -> 132 => (-1, 1, 0, 0)
# 241 -> 124 => (-1, 0, 0, 1)
# 241 -> 142 => (-1, 1, 0, 0)
# 241 -> 231 => (0, 0, -1, 1)

# julia> print_curve_classes(S)
# 124 -> 123: (-1, 1), Chern number: 0
# 132 -> 123: (1, 0), Chern number: 3
# 142 -> 124: (1, 0), Chern number: 3
# 142 -> 132: (0, 1), Chern number: 3
# 231 -> 123: (1, 0), Chern number: 3
# 231 -> 132: (1, 0), Chern number: 3
# 241 -> 124: (1, 0), Chern number: 3
# 241 -> 142: (1, 0), Chern number: 3
# 241 -> 231: (0, 1), Chern number: 3

R = root_system(:A, 3)
S = simple_roots(R)
G = generalized_gkm_schubert(R, S[1:1], "s1*s3*s2").self

println(G)
print_curve_classes(G)

beta = curve_class(G, "s3", "id") # Chern number 0 generator
gamma = curve_class(G, "s2", "id") # Chern number 2 generator

S1 = gkm_subgraph_from_vertices(G, ["id", "s2", "s3", "s3*s2"])
S2 = gkm_subgraph_from_vertices(G, ["id", "s1*s2", "s3", "s1*s3*s2"])
for d in 1:3
	println("[S1] * [S2] in q^($d*beta)")
	println(quantum_product(G, d*beta, poincare_dual(S1), poincare_dual(S2)))
end

# julia> integrate(poincare_dual(S2), edgeFromLabels(G, "id", "s3"))
# -1

# julia> integrate(poincare_dual(S1), edgeFromLabels(G, "id", "s3"))
# -1

# [S1] * [S2] in q^(1*beta)
# (t1*t2 - t1*t3 - t2*t3 + t3^2, 0, 0, 0, 0, t1*t2 - t1*t4 - t2*t4 + t4^2)
# [S1] * [S2] in q^(2*beta)
# (t1*t2 - t1*t3 - t2*t3 + t3^2, 0, 0, 0, 0, t1*t2 - t1*t4 - t2*t4 + t4^2)
# [S1] * [S2] in q^(3*beta)
# (t1*t2 - t1*t3 - t2*t3 + t3^2, 0, 0, 0, 0, t1*t2 - t1*t4 - t2*t4 + t4^2)

# julia> poincare_dual(gkm_subgraph_from_vertices(G, ["id", "s3"]))
# (t1*t2 - t1*t3 - t2*t3 + t3^2)*e[1] + (t1*t2 - t1*t4 - t2*t4 + t4^2)*e[6]

println()
println("Point quantum products:")

for k in 0:3
  for d in 0:3
 		println()
  	println("[id] * [id] in q^($d*beta + $k*gamma)")
  	println(quantum_product(G, d*beta+k*gamma, point_class(G, "id"), point_class(G, "id")))
  end
end

# julia> include("test/GW/examples/schubert/A3_a1_s1s3s2.jl")
# s2 -> id: (0, 1), Chern number: 3
# s1*s2 -> id: (0, 1), Chern number: 3
# s1*s2 -> s2: (0, 1), Chern number: 3
# s3*s2 -> s2: (1, 1), Chern number: 3
# s1*s3*s2 -> s1*s2: (1, 1), Chern number: 3
# s1*s3*s2 -> s3*s2: (0, 1), Chern number: 3
# s3 -> id: (1, 0), Chern number: 0
# s3 -> s3*s2: (0, 1), Chern number: 3
# s3 -> s1*s3*s2: (0, 1), Chern number: 3

# [id] * [id] in q^(0*beta + 0*gamma)
# (t1^2*t2^2*t3^2 - 2*t1^2*t2^2*t3*t4 + t1^2*t2^2*t4^2 - 2*t1^2*t2*t3^3 + 4*t1^2*t2*t3^2*t4 - 2*t1^2*t2*t3*t4^2 + t1^2*t3^4 - 2*t1^2*t3^3*t4 + t1^2*t3^2*t4^2 - 2*t1*t2^2*t3^3 + 4*t1*t2^2*t3^2*t4 - 2*t1*t2^2*t3*t4^2 + 4*t1*t2*t3^4 - 8*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - 2*t1*t3^5 + 4*t1*t3^4*t4 - 2*t1*t3^3*t4^2 + t2^2*t3^4 - 2*t2^2*t3^3*t4 + t2^2*t3^2*t4^2 - 2*t2*t3^5 + 4*t2*t3^4*t4 - 2*t2*t3^3*t4^2 + t3^6 - 2*t3^5*t4 + t3^4*t4^2, 0, 0, 0, 0, 0)

# [id] * [id] in q^(1*beta + 0*gamma)
# (t1^3*t2^3 - 3*t1^3*t2^2*t3 + 3*t1^3*t2*t3^2 - t1^3*t3^3 - 3*t1^2*t2^3*t3 + 9*t1^2*t2^2*t3^2 - 9*t1^2*t2*t3^3 + 3*t1^2*t3^4 + 3*t1*t2^3*t3^2 - 9*t1*t2^2*t3^3 + 9*t1*t2*t3^4 - 3*t1*t3^5 - t2^3*t3^3 + 3*t2^2*t3^4 - 3*t2*t3^5 + t3^6, 0, 0, 0, 0, t1^3*t2^3 - 2*t1^3*t2^2*t3 - t1^3*t2^2*t4 + t1^3*t2*t3^2 + 2*t1^3*t2*t3*t4 - t1^3*t3^2*t4 - 2*t1^2*t2^3*t3 - t1^2*t2^3*t4 + 4*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 + t1^2*t2^2*t4^2 - 2*t1^2*t2*t3^3 - 5*t1^2*t2*t3^2*t4 - 2*t1^2*t2*t3*t4^2 + 2*t1^2*t3^3*t4 + t1^2*t3^2*t4^2 + t1*t2^3*t3^2 + 2*t1*t2^3*t3*t4 - 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 - 2*t1*t2^2*t3*t4^2 + t1*t2*t3^4 + 4*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - t1*t3^4*t4 - 2*t1*t3^3*t4^2 - t2^3*t3^2*t4 + 2*t2^2*t3^3*t4 + t2^2*t3^2*t4^2 - t2*t3^4*t4 - 2*t2*t3^3*t4^2 + t3^4*t4^2)

# [id] * [id] in q^(2*beta + 0*gamma)
# (t1^3*t2^3 - 3*t1^3*t2^2*t3 + 3*t1^3*t2*t3^2 - t1^3*t3^3 - 3*t1^2*t2^3*t3 + 9*t1^2*t2^2*t3^2 - 9*t1^2*t2*t3^3 + 3*t1^2*t3^4 + 3*t1*t2^3*t3^2 - 9*t1*t2^2*t3^3 + 9*t1*t2*t3^4 - 3*t1*t3^5 - t2^3*t3^3 + 3*t2^2*t3^4 - 3*t2*t3^5 + t3^6, 0, 0, 0, 0, t1^3*t2^3 - 2*t1^3*t2^2*t3 - t1^3*t2^2*t4 + t1^3*t2*t3^2 + 2*t1^3*t2*t3*t4 - t1^3*t3^2*t4 - 2*t1^2*t2^3*t3 - t1^2*t2^3*t4 + 4*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 + t1^2*t2^2*t4^2 - 2*t1^2*t2*t3^3 - 5*t1^2*t2*t3^2*t4 - 2*t1^2*t2*t3*t4^2 + 2*t1^2*t3^3*t4 + t1^2*t3^2*t4^2 + t1*t2^3*t3^2 + 2*t1*t2^3*t3*t4 - 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 - 2*t1*t2^2*t3*t4^2 + t1*t2*t3^4 + 4*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - t1*t3^4*t4 - 2*t1*t3^3*t4^2 - t2^3*t3^2*t4 + 2*t2^2*t3^3*t4 + t2^2*t3^2*t4^2 - t2*t3^4*t4 - 2*t2*t3^3*t4^2 + t3^4*t4^2)

# [id] * [id] in q^(3*beta + 0*gamma)
# (t1^3*t2^3 - 3*t1^3*t2^2*t3 + 3*t1^3*t2*t3^2 - t1^3*t3^3 - 3*t1^2*t2^3*t3 + 9*t1^2*t2^2*t3^2 - 9*t1^2*t2*t3^3 + 3*t1^2*t3^4 + 3*t1*t2^3*t3^2 - 9*t1*t2^2*t3^3 + 9*t1*t2*t3^4 - 3*t1*t3^5 - t2^3*t3^3 + 3*t2^2*t3^4 - 3*t2*t3^5 + t3^6, 0, 0, 0, 0, t1^3*t2^3 - 2*t1^3*t2^2*t3 - t1^3*t2^2*t4 + t1^3*t2*t3^2 + 2*t1^3*t2*t3*t4 - t1^3*t3^2*t4 - 2*t1^2*t2^3*t3 - t1^2*t2^3*t4 + 4*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 + t1^2*t2^2*t4^2 - 2*t1^2*t2*t3^3 - 5*t1^2*t2*t3^2*t4 - 2*t1^2*t2*t3*t4^2 + 2*t1^2*t3^3*t4 + t1^2*t3^2*t4^2 + t1*t2^3*t3^2 + 2*t1*t2^3*t3*t4 - 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 - 2*t1*t2^2*t3*t4^2 + t1*t2*t3^4 + 4*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - t1*t3^4*t4 - 2*t1*t3^3*t4^2 - t2^3*t3^2*t4 + 2*t2^2*t3^3*t4 + t2^2*t3^2*t4^2 - t2*t3^4*t4 - 2*t2*t3^3*t4^2 + t3^4*t4^2)

# [id] * [id] in q^(0*beta + 1*gamma)
# (t1*t3^2 - 2*t1*t3*t4 + t1*t4^2 + t2*t3^2 - 2*t2*t3*t4 + t2*t4^2 - 2*t3^3 + 4*t3^2*t4 - 2*t3*t4^2, t1*t3^2 - 2*t1*t3*t4 + t1*t4^2 - t3^3 + 2*t3^2*t4 - t3*t4^2, t2*t3^2 - 2*t2*t3*t4 + t2*t4^2 - t3^3 + 2*t3^2*t4 - t3*t4^2, 0, 0, 0)

# [id] * [id] in q^(1*beta + 1*gamma)
# (2*t1^2*t2 - t1^2*t3 - t1^2*t4 + 2*t1*t2^2 - 7*t1*t2*t3 - t1*t2*t4 + 3*t1*t3^2 + 3*t1*t3*t4 - t2^2*t3 - t2^2*t4 + 3*t2*t3^2 + 3*t2*t3*t4 - t3^3 - 3*t3^2*t4, t1^2*t2 - t1^2*t4 - 2*t1*t2*t3 + 2*t1*t3*t4 + t2*t3^2 - t3^2*t4, t1*t2^2 - 2*t1*t2*t3 + t1*t3^2 - t2^2*t4 + 2*t2*t3*t4 - t3^2*t4, t1^2*t2 - t1^2*t3 - 2*t1*t2*t3 + 2*t1*t3^2 + t2*t3^2 - t3^3, t1*t2^2 - 2*t1*t2*t3 + t1*t3^2 - t2^2*t3 + 2*t2*t3^2 - t3^3, 2*t1^2*t2 - t1^2*t3 - t1^2*t4 + 2*t1*t2^2 - 5*t1*t2*t3 - 3*t1*t2*t4 + 2*t1*t3^2 + 3*t1*t3*t4 + t1*t4^2 - t2^2*t3 - t2^2*t4 + 2*t2*t3^2 + 3*t2*t3*t4 + t2*t4^2 - t3^3 - t3^2*t4 - 2*t3*t4^2)

# [id] * [id] in q^(2*beta + 1*gamma)
# (0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 1*gamma)
# (0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(0*beta + 2*gamma)
# (0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(1*beta + 2*gamma)
# (1, 1, 1, 1, 1, 1)

# [id] * [id] in q^(2*beta + 2*gamma)
# (0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 2*gamma)
# (0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(0*beta + 3*gamma)
# (0, 0, 0, 0, 0, 0)