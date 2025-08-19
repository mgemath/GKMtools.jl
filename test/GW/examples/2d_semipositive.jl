# In this file, we calculate certain Gromov--Witten invariants and quantum products of
# the Hirzebruch 2-surface \Sigma_2.
# It is an interesting example because the Chern classes of its weights are non-negative with
# a unique edge of Chern class zero, yielding infinitely many contributions to the quantum product.

# Hirzebruch 2:
G = gkm_2d([1 0; 1 1; -1 0; 1 -1])

# julia> print_curve_classes(G)
# 2 -> 1: (-2, 1), Chern number: 0
# 3 -> 2: (1, 0), Chern number: 2
# 4 -> 1: (1, 0), Chern number: 2
# 4 -> 3: (0, 1), Chern number: 4

b1 = curve_class(G, Edge(2, 3))
b0 = curve_class(G, Edge(1, 2))

for d in 1:5
  println("Class $d*b0 = $(d*b0), 0 points -> $(gromov_witten(G, d*b0, 0, class_one()))")
end

# Class 1*b0 = (-2, 1), 0 points -> t2
# Class 2*b0 = (-4, 2), 0 points -> 1//8*t2
# Class 3*b0 = (-6, 3), 0 points -> 1//27*t2
# Class 4*b0 = (-8, 4), 0 points -> 1//64*t2
# Class 5*b0 = (-10, 5), 0 points -> 1//125*t2

# julia> for d in 1:5
#        println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 1, ev(1, point_class(1, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> -t1*t2 + t2^2
# d=2, d*b0 = (-4, 2) -> -1//4*t1*t2 + 1//4*t2^2
# d=3, d*b0 = (-6, 3) -> -1//9*t1*t2 + 1//9*t2^2
# d=4, d*b0 = (-8, 4) -> -1//16*t1*t2 + 1//16*t2^2
# d=5, d*b0 = (-10, 5) -> -1//25*t1*t2 + 1//25*t2^2

# julia> for d in 1:5
#          println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 2, ev(1, point_class(1, G)) * ev(2, point_class(2, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> -t1^2*t2 + t2^3
# d=2, d*b0 = (-4, 2) -> -1//2*t1^2*t2 + 1//2*t2^3
# d=3, d*b0 = (-6, 3) -> -1//3*t1^2*t2 + 1//3*t2^3
# d=4, d*b0 = (-8, 4) -> -1//4*t1^2*t2 + 1//4*t2^3
# d=5, d*b0 = (-10, 5) -> -1//5*t1^2*t2 + 1//5*t2^3

# julia> for d in 1:5
#        println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 3, ev(1, point_class(1, G)) * ev(2, point_class(1, G)) * ev(3, point_class(2, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> t1^3*t2 - t1^2*t2^2 - t1*t2^3 + t2^4
# d=2, d*b0 = (-4, 2) -> t1^3*t2 - t1^2*t2^2 - t1*t2^3 + t2^4
# d=3, d*b0 = (-6, 3) -> t1^3*t2 - t1^2*t2^2 - t1*t2^3 + t2^4
# d=4, d*b0 = (-8, 4) -> t1^3*t2 - t1^2*t2^2 - t1*t2^3 + t2^4

# julia> for d in 1:5
#        println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 4, ev(1, point_class(1, G)) * ev(2, point_class(1, G)) * ev(3, point_class(2, G)) * ev(4, point_class(2, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> t1^4*t2 - 2*t1^2*t2^3 + t2^5
# d=2, d*b0 = (-4, 2) -> 2*t1^4*t2 - 4*t1^2*t2^3 + 2*t2^5
# d=3, d*b0 = (-6, 3) -> 3*t1^4*t2 - 6*t1^2*t2^3 + 3*t2^5
# d=4, d*b0 = (-8, 4) -> 4*t1^4*t2 - 8*t1^2*t2^3 + 4*t2^5

# julia> for d in 1:4
#        println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 5, ev(1, point_class(1, G)) * ev(2, point_class(1, G)) * ev(3, point_class(2, G)) * ev(4, point_class(2, G)) * ev(5, point_class(1, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> -t1^5*t2 + t1^4*t2^2 + 2*t1^3*t2^3 - 2*t1^2*t2^4 - t1*t2^5 + t2^6
# d=2, d*b0 = (-4, 2) -> -4*t1^5*t2 + 4*t1^4*t2^2 + 8*t1^3*t2^3 - 8*t1^2*t2^4 - 4*t1*t2^5 + 4*t2^6
# d=3, d*b0 = (-6, 3) -> -9*t1^5*t2 + 9*t1^4*t2^2 + 18*t1^3*t2^3 - 18*t1^2*t2^4 - 9*t1*t2^5 + 9*t2^6
# d=4, d*b0 = (-8, 4) -> -16*t1^5*t2 + 16*t1^4*t2^2 + 32*t1^3*t2^3 - 32*t1^2*t2^4 - 16*t1*t2^5 + 16*t2^6


n = valency(G)
maxChernNumber = 2*n
nv = n_vertices(G.g)
evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
println("Starting integration:")

for i in 0:5
  for j in 0:5
    i == 0 && j == 0 && continue
    beta = i*b1 + j*b0
    println("Calculating structure constants for $beta, Chern number $(chern_number(G, beta)):")
    QH_structure_constants(G, beta; P_input = P_input, show_progress = true)
  end
end

# Zero chern class contributions seem to be constant!

# julia> for i in 1:5
# println("S[$i * b0] == S[b0] -> $(S[i*b0]==S[b0])")
# end
# S[1 * b0] == S[b0] -> true
# S[2 * b0] == S[b0] -> true
# S[3 * b0] == S[b0] -> true
# S[4 * b0] == S[b0] -> true
# S[5 * b0] == S[b0] -> true