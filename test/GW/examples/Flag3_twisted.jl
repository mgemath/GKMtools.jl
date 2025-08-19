# This file experiments with the twisted flag variety of Tolman/Eschenburg/Woodward.
# Again, it is interesting that it contains a unique edge with zero chern number, while all others are positive.

G = gkm_3d_twisted_flag()

# Curve classes:
# julia> print_curve_classes(ans)
# 2 -> 1: (0, 1), Chern number: 4
# 3 -> 2: (-1, 1), Chern number: 2
# 4 -> 1: (1, 0), Chern number: 2
# 4 -> 3: (-2, 1), Chern number: 0
# 5 -> 2: (1, 0), Chern number: 2
# 5 -> 4: (-1, 1), Chern number: 2
# 6 -> 1: (1, 1), Chern number: 6
# 6 -> 3: (1, 0), Chern number: 2
# 6 -> 5: (0, 1), Chern number: 4

b1 = curve_class(G, Edge(2, 3)) # (-1, 1), c1 = 2
b2 = curve_class(G, Edge(1, 4)) # (1, 0), c1 = 2
b0 = curve_class(G, Edge(3, 4)) # (-2, 1), c1 = 0

for d in 1:5
  println("Class $d*b0 = $(d*b0), 0 points -> $(gromov_witten(G, d*b0, 0, class_one(); show_bar=true))")
end

# Class 1*b0 = (-2, 1), 0 points -> 1
# Class 2*b0 = (-4, 2), 0 points -> -7//8
# Class 3*b0 = (-6, 3), 0 points -> 55//27
# Class 4*b0 = (-8, 4), 0 points -> -455//64
# Class 5*b0 = (-10, 5), 0 points -> 3876//

# julia> for d in 1:6
#        println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 1, ev(1, point_class(3, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> -2*t1^2 + t1*t2
# d=2, d*b0 = (-4, 2) -> 7//2*t1^2 - 7//4*t1*t2
# d=3, d*b0 = (-6, 3) -> -110//9*t1^2 + 55//9*t1*t2
# d=4, d*b0 = (-8, 4) -> 455//8*t1^2 - 455//16*t1*t2
# d=5, d*b0 = (-10, 5) -> -7752//25*t1^2 + 3876//25*t1*t2

# julia> for d in 1:5
#          println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 2, ev(1, point_class(3, G)) * ev(2, point_class(3, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> 4*t1^4 - 4*t1^3*t2 + t1^2*t2^2
# d=2, d*b0 = (-4, 2) -> -14*t1^4 + 14*t1^3*t2 - 7//2*t1^2*t2^2
# d=3, d*b0 = (-6, 3) -> 220//3*t1^4 - 220//3*t1^3*t2 + 55//3*t1^2*t2^2
# d=4, d*b0 = (-8, 4) -> -455*t1^4 + 455*t1^3*t2 - 455//4*t1^2*t2^2
# d=5, d*b0 = (-10, 5) -> 15504//5*t1^4 - 15504//5*t1^3*t2 + 3876//5*t1^2*t2^2

# julia> for d in 1:5
#          println("d=$d, d*b0 = $(d*b0) -> $(gromov_witten(G, d*b0, 4, ev(1, point_class(3, G)) * ev(2, point_class(3, G)) * ev(3, point_class(3, G)) * ev(4, point_class(3, G))))")
#        end
# d=1, d*b0 = (-2, 1) -> 16*t1^8 - 32*t1^7*t2 + 24*t1^6*t2^2 - 8*t1^5*t2^3 + t1^4*t2^4
# d=2, d*b0 = (-4, 2) -> -224*t1^8 + 448*t1^7*t2 - 336*t1^6*t2^2 + 112*t1^5*t2^3 - 14*t1^4*t2^4
# d=3, d*b0 = (-6, 3) -> 2640*t1^8 - 5280*t1^7*t2 + 3960*t1^6*t2^2 - 1320*t1^5*t2^3 + 165*t1^4*t2^4
# d=4, d*b0 = (-8, 4) -> -29120*t1^8 + 58240*t1^7*t2 - 43680*t1^6*t2^2 + 14560*t1^5*t2^3 - 1820*t1^4*t2^4
# d=5, d*b0 = (-10, 5) -> 310080*t1^8 - 620160*t1^7*t2 + 465120*t1^6*t2^2 - 155040*t1^5*t2^3 + 19380*t1^4*t2^4

n = valency(G)
maxChernNumber = 2*n
nv = n_vertices(G.g)
evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
println("Starting integration:")

betas = [zero(parent(b1)), b1, b2, 2*b1, b1+b2, 2*b2, 3*b1, 2*b1+b2, b1+2*b2, 3*b2]

for b in betas
  continue
  for i in 0:6
    beta = b + i*b0
    println("Calculating structure constants for $beta, Chern number $(chern_number(G, beta)):")
    QH_structure_constants(G, beta; P_input = P_input, show_progress = true)
  end
end

for i in 0:5
  continue
  beta = 2*b2 + i*b0
  println("Calculating structure constants for $beta, Chern number $(chern_number(G, beta)):")
  QH_structure_constants(G, beta; P_input = P_input, show_progress = true)
end

single_P_input = evClasses[1, 1]*evClasses[2, 1]*evClasses[3, 1]

for i in 5:6
  beta = 3*b2 + i*b0
  println("Calculating structure constants for $beta, Chern number $(chern_number(G, beta)):")
  println("Term is: $(gromov_witten(G, beta, 3, single_P_input; show_bar=true))")
end

function check_flag3_phenomenon(S)
  m = [1, -7, 55, -455, 3876, -33649]
  for i in 1:6
    t = S[i*b0][3:4,3:4,3:4] == m[i] .* S[b0][3:4,3:4,3:4]
    println("Test for i=$i passes -> $t")
  end
end

#TODO: iterate through multiplicities of the c1=0 curve class.
#Interesting Q: After specialising to t=0, does QH still see the difference between F3 and F3_twisted?
