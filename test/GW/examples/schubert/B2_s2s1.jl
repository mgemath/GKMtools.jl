# This file deals with the smooth Schubert variety in G/B in type B2, given by s1*s2.

R = root_system(:B, 2)
S = simple_roots(R)
G = generalized_gkm_schubert(R, S[1:0], "s1*s2").self

print_curve_classes(G)

beta = curve_class(G, "s2", "id") # Chern number 0 generator
gamma = curve_class(G, "id", "s1") # Chern number 2 generator

for k in 0:3
	for d in 0:3
		println()
		println("[s1] * [s1] in q^($d*beta + $k*gamma)")
		println(quantum_product(G, d*beta+k*gamma, point_class(G, "id"), point_class(G, "id")))
	end
end

# julia> include("test/GW/examples/schubert/B2_s2s1.jl")
# s1 -> id: (-1, 1), Chern number: 1
# s2 -> id: (1, 0), Chern number: 2
# s1*s2 -> s1: (1, 0), Chern number: 2
# s1*s2 -> s2: (0, 1), Chern number: 3

# [s1] * [s1] in q^(0*beta + 0*gamma)
# (t1^2*t2^2 - 2*t1*t2^3 + t2^4, 0, 0, 0)

# [s1] * [s1] in q^(1*beta + 0*gamma)
# (t1^2 - 2*t1*t2 + t2^2, 0, t1^2 - 2*t1*t2 + t2^2, 0)

# [s1] * [s1] in q^(2*beta + 0*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*beta + 0*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*beta + 1*gamma)
# (t2^3, t1*t2^2, 0, 0)

# [s1] * [s1] in q^(1*beta + 1*gamma)
# (t1 + t2, t1 + t2, t1, t2)

# [s1] * [s1] in q^(2*beta + 1*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*beta + 1*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*beta + 2*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(1*beta + 2*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(2*beta + 2*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*beta + 2*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*beta + 3*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(1*beta + 3*gamma)
# (0, 0, 0, 0)

# [s1] * [s1] in q^(2*beta + 3*gamma)
# (0, 0, 0, 0)
