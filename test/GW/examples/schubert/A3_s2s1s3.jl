# In this file we study the smooth Schubert variety in type A3 in G/B given by s2*s1*s3.

R = root_system(:A, 3)
S = simple_roots(R)
G = generalized_gkm_schubert(R, S[1:0], "s2*s1*s3").self

print_curve_classes(G)

beta = curve_class(G, "s2", "id") # Chern number 0 generator
gamma2 = curve_class(G, "s1", "id") # Chern number 2 generator
gamma3 = curve_class(G, "id", "s3") # Chern number 2 generator

S1 = gkm_subgraph_from_vertices(G, ["id", "s2", "s1", "s2*s1"])
S2 = gkm_subgraph_from_vertices(G, ["id", "s2", "s3", "s2*s3"])

# julia> integrate(poincare_dual(S1), G, edgeFromLabels(G, "id", "s2"))
# -1

# julia> integrate(poincare_dual(S2), G, edgeFromLabels(G, "id", "s2"))
# -1

for k2 in 0:3
  for k3 in 0:3
	  for d in 0:3
	   println()
	   println("[id] * [id] in q^($d*beta + $k2*gamma2 + $k3*gamma3)")
	   println(quantum_product(G, d*beta+k2*gamma2+k3*gamma3, point_class(G, "id"), point_class(G, "id")))
    end
  end
end

# julia> include("test/GW/examples/schubert/A3_s2s1s3.jl")
# s1 -> id: (0, 0, 1), Chern number: 2
# s2*s1 -> s1: (-1, 1, 1), Chern number: 2
# s2 -> id: (-1, 1, 0), Chern number: 0
# s2 -> s2*s1: (0, 0, 1), Chern number: 2
# s3 -> id: (1, 0, 0), Chern number: 2
# s1*s3 -> s1: (1, 0, 0), Chern number: 2
# s1*s3 -> s3: (0, 0, 1), Chern number: 2
# s2*s1*s3 -> s2*s1: (1, 0, 0), Chern number: 2
# s2*s1*s3 -> s1*s3: (0, 1, 1), Chern number: 4
# s2*s3 -> s2: (1, 0, 0), Chern number: 2
# s2*s3 -> s3: (0, 1, 0), Chern number: 2
# s2*s3 -> s2*s1*s3: (0, 0, 1), Chern number: 2

# [id] * [id] in q^(0*beta + 0*gamma2 + 0*gamma3)
# (t1^2*t2^2*t3^2 - 2*t1^2*t2^2*t3*t4 + t1^2*t2^2*t4^2 - 2*t1^2*t2*t3^3 + 4*t1^2*t2*t3^2*t4 - 2*t1^2*t2*t3*t4^2 + t1^2*t3^4 - 2*t1^2*t3^3*t4 + t1^2*t3^2*t4^2 - 2*t1*t2^3*t3^2 + 4*t1*t2^3*t3*t4 - 2*t1*t2^3*t4^2 + 4*t1*t2^2*t3^3 - 8*t1*t2^2*t3^2*t4 + 4*t1*t2^2*t3*t4^2 - 2*t1*t2*t3^4 + 4*t1*t2*t3^3*t4 - 2*t1*t2*t3^2*t4^2 + t2^4*t3^2 - 2*t2^4*t3*t4 + t2^4*t4^2 - 2*t2^3*t3^3 + 4*t2^3*t3^2*t4 - 2*t2^3*t3*t4^2 + t2^2*t3^4 - 2*t2^2*t3^3*t4 + t2^2*t3^2*t4^2, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(1*beta + 0*gamma2 + 0*gamma3)
# (t1^3*t3^3 - 3*t1^3*t3^2*t4 + 3*t1^3*t3*t4^2 - t1^3*t4^3 - 3*t1^2*t2*t3^3 + 9*t1^2*t2*t3^2*t4 - 9*t1^2*t2*t3*t4^2 + 3*t1^2*t2*t4^3 + 3*t1*t2^2*t3^3 - 9*t1*t2^2*t3^2*t4 + 9*t1*t2^2*t3*t4^2 - 3*t1*t2^2*t4^3 - t2^3*t3^3 + 3*t2^3*t3^2*t4 - 3*t2^3*t3*t4^2 + t2^3*t4^3, 0, 0, t1^3*t2*t3^2 - 2*t1^3*t2*t3*t4 + t1^3*t2*t4^2 - t1^3*t3^2*t4 + 2*t1^3*t3*t4^2 - t1^3*t4^3 - 2*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 - 2*t1^2*t2^2*t4^2 - t1^2*t2*t3^3 + 4*t1^2*t2*t3^2*t4 - 5*t1^2*t2*t3*t4^2 + 2*t1^2*t2*t4^3 + t1^2*t3^3*t4 - 2*t1^2*t3^2*t4^2 + t1^2*t3*t4^3 + t1*t2^3*t3^2 - 2*t1*t2^3*t3*t4 + t1*t2^3*t4^2 + 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 + 4*t1*t2^2*t3*t4^2 - t1*t2^2*t4^3 - 2*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - 2*t1*t2*t3*t4^3 - t2^3*t3^3 + 2*t2^3*t3^2*t4 - t2^3*t3*t4^2 + t2^2*t3^3*t4 - 2*t2^2*t3^2*t4^2 + t2^2*t3*t4^3, 0, 0, 0, 0)

# [id] * [id] in q^(2*beta + 0*gamma2 + 0*gamma3)
# (t1^3*t3^3 - 3*t1^3*t3^2*t4 + 3*t1^3*t3*t4^2 - t1^3*t4^3 - 3*t1^2*t2*t3^3 + 9*t1^2*t2*t3^2*t4 - 9*t1^2*t2*t3*t4^2 + 3*t1^2*t2*t4^3 + 3*t1*t2^2*t3^3 - 9*t1*t2^2*t3^2*t4 + 9*t1*t2^2*t3*t4^2 - 3*t1*t2^2*t4^3 - t2^3*t3^3 + 3*t2^3*t3^2*t4 - 3*t2^3*t3*t4^2 + t2^3*t4^3, 0, 0, t1^3*t2*t3^2 - 2*t1^3*t2*t3*t4 + t1^3*t2*t4^2 - t1^3*t3^2*t4 + 2*t1^3*t3*t4^2 - t1^3*t4^3 - 2*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 - 2*t1^2*t2^2*t4^2 - t1^2*t2*t3^3 + 4*t1^2*t2*t3^2*t4 - 5*t1^2*t2*t3*t4^2 + 2*t1^2*t2*t4^3 + t1^2*t3^3*t4 - 2*t1^2*t3^2*t4^2 + t1^2*t3*t4^3 + t1*t2^3*t3^2 - 2*t1*t2^3*t3*t4 + t1*t2^3*t4^2 + 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 + 4*t1*t2^2*t3*t4^2 - t1*t2^2*t4^3 - 2*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - 2*t1*t2*t3*t4^3 - t2^3*t3^3 + 2*t2^3*t3^2*t4 - t2^3*t3*t4^2 + t2^2*t3^3*t4 - 2*t2^2*t3^2*t4^2 + t2^2*t3*t4^3, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 0*gamma2 + 0*gamma3)
# (t1^3*t3^3 - 3*t1^3*t3^2*t4 + 3*t1^3*t3*t4^2 - t1^3*t4^3 - 3*t1^2*t2*t3^3 + 9*t1^2*t2*t3^2*t4 - 9*t1^2*t2*t3*t4^2 + 3*t1^2*t2*t4^3 + 3*t1*t2^2*t3^3 - 9*t1*t2^2*t3^2*t4 + 9*t1*t2^2*t3*t4^2 - 3*t1*t2^2*t4^3 - t2^3*t3^3 + 3*t2^3*t3^2*t4 - 3*t2^3*t3*t4^2 + t2^3*t4^3, 0, 0, t1^3*t2*t3^2 - 2*t1^3*t2*t3*t4 + t1^3*t2*t4^2 - t1^3*t3^2*t4 + 2*t1^3*t3*t4^2 - t1^3*t4^3 - 2*t1^2*t2^2*t3^2 + 4*t1^2*t2^2*t3*t4 - 2*t1^2*t2^2*t4^2 - t1^2*t2*t3^3 + 4*t1^2*t2*t3^2*t4 - 5*t1^2*t2*t3*t4^2 + 2*t1^2*t2*t4^3 + t1^2*t3^3*t4 - 2*t1^2*t3^2*t4^2 + t1^2*t3*t4^3 + t1*t2^3*t3^2 - 2*t1*t2^3*t3*t4 + t1*t2^3*t4^2 + 2*t1*t2^2*t3^3 - 5*t1*t2^2*t3^2*t4 + 4*t1*t2^2*t3*t4^2 - t1*t2^2*t4^3 - 2*t1*t2*t3^3*t4 + 4*t1*t2*t3^2*t4^2 - 2*t1*t2*t3*t4^3 - t2^3*t3^3 + 2*t2^3*t3^2*t4 - t2^3*t3*t4^2 + t2^2*t3^3*t4 - 2*t2^2*t3^2*t4^2 + t2^2*t3*t4^3, 0, 0, 0, 0)

# [id] * [id] in q^(0*beta + 0*gamma2 + 1*gamma3)
# (t1^2*t2^2 - 2*t1^2*t2*t3 + t1^2*t3^2 - 2*t1*t2^3 + 4*t1*t2^2*t3 - 2*t1*t2*t3^2 + t2^4 - 2*t2^3*t3 + t2^2*t3^2, 0, 0, 0, t1^2*t2^2 - 2*t1^2*t2*t3 + t1^2*t3^2 - 2*t1*t2^3 + 4*t1*t2^2*t3 - 2*t1*t2*t3^2 + t2^4 - 2*t2^3*t3 + t2^2*t3^2, 0, 0, 0)

# [id] * [id] in q^(1*beta + 0*gamma2 + 1*gamma3)
# (t1^3*t2 + t1^3*t3 - 2*t1^3*t4 - 3*t1^2*t2^2 - 3*t1^2*t2*t3 + 6*t1^2*t2*t4 + 3*t1*t2^3 + 3*t1*t2^2*t3 - 6*t1*t2^2*t4 - t2^4 - t2^3*t3 + 2*t2^3*t4, 0, 0, t1^3*t2 + t1^3*t3 - 2*t1^3*t4 - 2*t1^2*t2^2 - 3*t1^2*t2*t3 + 4*t1^2*t2*t4 - t1^2*t3^2 + 2*t1^2*t3*t4 + t1*t2^3 + 3*t1*t2^2*t3 - 2*t1*t2^2*t4 + 2*t1*t2*t3^2 - 4*t1*t2*t3*t4 - t2^3*t3 - t2^2*t3^2 + 2*t2^2*t3*t4, t1^3*t2 - t1^3*t4 - 3*t1^2*t2^2 + 3*t1^2*t2*t4 + 3*t1*t2^3 - 3*t1*t2^2*t4 - t2^4 + t2^3*t4, 0, 0, t1^3*t3 - t1^3*t4 - 2*t1^2*t2*t3 + 2*t1^2*t2*t4 - t1^2*t3^2 + t1^2*t3*t4 + t1*t2^2*t3 - t1*t2^2*t4 + 2*t1*t2*t3^2 - 2*t1*t2*t3*t4 - t2^2*t3^2 + t2^2*t3*t4)

# [id] * [id] in q^(2*beta + 0*gamma2 + 1*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 0*gamma2 + 1*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(0*beta + 0*gamma2 + 2*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(1*beta + 0*gamma2 + 2*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(2*beta + 0*gamma2 + 2*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 0*gamma2 + 2*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(0*beta + 0*gamma2 + 3*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(1*beta + 0*gamma2 + 3*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(2*beta + 0*gamma2 + 3*gamma3)
# (0, 0, 0, 0, 0, 0, 0, 0)

# [id] * [id] in q^(3*beta + 0*gamma2 + 3*gamma3)
# Killed