# This file deals with the smooth Schubert variety in G/B in type G2, given by s1*s2.

R = root_system(:G, 2)
S = simple_roots(R)
G = generalized_gkm_schubert(R, S[1:0], "s1*s2").self

print_curve_classes(G)

alpha = curve_class(G, "s1", "id") # Chern number -1 generator
gamma = curve_class(G, "id", "s2") # Chern number 2 generator

for k in 0:3
	for d in 0:3
		println()
		println("[s1] * [s1] in q^($d*alpha + $k*gamma), chern number: $(-d+2*k)")
		println(quantum_product(G, d*beta+k*gamma, point_class(G, "id"), point_class(G, "id")))
	end
end

# julia> include("test/GW/examples/schubert/G2_s1s2.jl")
# s1 -> id: (-3, 1), Chern number: -1
# s2 -> id: (1, 0), Chern number: 2
# s1*s2 -> s1: (1, 0), Chern number: 2
# s1*s2 -> s2: (0, 1), Chern number: 5

# [s1] * [s1] in q^(0*alpha + 0*gamma), chern number: 0
# (4*t1^4 - 12*t1^3*t2 - 4*t1^3*t3 + 13*t1^2*t2^2 + 10*t1^2*t2*t3 + t1^2*t3^2 - 6*t1*t2^3 - 8*t1*t2^2*t3 - 2*t1*t2*t3^2 + t2^4 + 2*t2^3*t3 + t2^2*t3^2, 0, 0, 0)

# [s1] * [s1] in q^(1*alpha + 0*gamma), chern number: -1
# (t1^2 - 2*t1*t2 + t2^2, 0, t1^2 - 2*t1*t2 + t2^2, 0)

# [s1] * [s1] in q^(2*alpha + 0*gamma), chern number: -2
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*alpha + 0*gamma), chern number: -3
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*alpha + 1*gamma), chern number: 2
# (t1^2 - 2*t1*t2 + t2^2, 0, t1^2 - 2*t1*t2 + t2^2, 0)

# [s1] * [s1] in q^(1*alpha + 1*gamma), chern number: 1
# (0, 0, 0, 0)

# [s1] * [s1] in q^(2*alpha + 1*gamma), chern number: 0
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*alpha + 1*gamma), chern number: -1
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*alpha + 2*gamma), chern number: 4
# (0, 0, 0, 0)

# [s1] * [s1] in q^(1*alpha + 2*gamma), chern number: 3
# (0, 0, 0, 0)

# [s1] * [s1] in q^(2*alpha + 2*gamma), chern number: 2
# (0, 0, 0, 0)

# [s1] * [s1] in q^(3*alpha + 2*gamma), chern number: 1
# (0, 0, 0, 0)

# [s1] * [s1] in q^(0*alpha + 3*gamma), chern number: 6
# (0, 0, 0, 0)

# [s1] * [s1] in q^(1*alpha + 3*gamma), chern number: 5
# (0, 0, 0, 0)

# [s1] * [s1] in q^(2*alpha + 3*gamma), chern number: 4
# (0, 0, 0, 0)