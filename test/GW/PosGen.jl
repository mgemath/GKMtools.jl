H = load_H()

# P2 = projective_space(GKM_graph, 2)
# beta = curve_class(P2, Edge(1, 2))
# 
# P2_g1_m3 = GKMtools.gromov_witten_pos_gen(P2, 1*beta, 3, 1, prod([ev(i, point_class(P2, 2)) for i in 1:3]), H)
# P2_g1_m6 = GKMtools.gromov_witten_pos_gen(P2, 2*beta, 6, 1, prod([ev(i, point_class(P2, 2)) for i in 1:6]), H)
# P2_g2_m4 = GKMtools.gromov_witten_pos_gen(P2, 1*beta, 4, 2, prod([ev(i, point_class(P2, 2)) for i in 1:4]), H)
# P2_g2_m7 = GKMtools.gromov_witten_pos_gen(P2, 2*beta, 7, 2, prod([ev(i, point_class(P2, 2)) for i in 1:7]), H)
# 
# println([P2_g1_m3, P2_g1_m6, P2_g2_m4, P2_g2_m7])

F = gkm_3d_twisted_flag()
gamma = curve_class(F, Edge(3 ,4))
for g in 0:3
  for d in 1:4
    gw = GKMtools.gromov_witten_pos_gen(F, d*gamma, 0, g, class_one(), H)[1]
    #println("Twisted flag, g=$g, d=$d: $(gw)")
    println("Twisted flag, g=$g, d=$d: denominator is $(factor(denominator(gw)))")
  end
end

#### Denominators WITHOUT experimental mode:
#### Exponent seems to be 2g + 2d - 2
#### (t1 - t2) is the edge weight.
# Twisted flag, g=0, d=1: denominator is 1
# Twisted flag, g=0, d=2: denominator is 1
# Twisted flag, g=0, d=3: denominator is 1
# Twisted flag, g=0, d=4: denominator is 1
# Twisted flag, g=1, d=1: denominator is 1
# Twisted flag, g=1, d=2: denominator is 1 * (t1 - t2)^4
# Twisted flag, g=1, d=3: denominator is 1 * (t1 - t2)^6
# Twisted flag, g=1, d=4: denominator is 1 * (t1 - t2)^8
# Twisted flag, g=2, d=1: denominator is 1
# Twisted flag, g=2, d=2: denominator is 1 * (t1 - t2)^6
# Twisted flag, g=2, d=3: denominator is 1 * (t1 - t2)^8
# Twisted flag, g=2, d=4: denominator is 1 * (t1 - t2)^10
# Twisted flag, g=3, d=1: denominator is 1
# Twisted flag, g=3, d=2: denominator is 1 * (t1 - t2)^8
# Twisted flag, g=3, d=3: denominator is 1 * (t1 - t2)^10
# Twisted flag, g=3, d=4: denominator is 1 * (t1 - t2)^12
####

#### Output WITH experimental mode (i.e. line 18 in Euler_pos_gen.jl):
#
# AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}[0, 0, 0, 0]
# Twisted flag, g=0, d=1: 1
# Twisted flag, g=0, d=2: -7//8
# Twisted flag, g=0, d=3: 55//27
# Twisted flag, g=0, d=4: -455//64
# Twisted flag, g=0, d=5: 3876//125
# Twisted flag, g=1, d=1: 1//12
# Twisted flag, g=1, d=2: (-1//24*t1^2 + 23//96*t1*t2 - 1//24*t2^2)//(t1^2 - 2*t1*t2 + t2^2)
# Twisted flag, g=1, d=3: (-557//2592*t1^6 - 499//2592*t1^5*t2 + 5101//2592*t1^4*t2^2 - 4117//1296*t1^3*t2^3 + 5101//2592*t1^2*t2^4 - 499//2592*t1*t2^5 - 557//2592*t2^6)//(t1^6 - 6*t1^5*t2 + 15*t1^4*t2^2 - 20*t1^3*t2^3 + 15*t1^2*t2^4 - 6*t1*t2^5 + t2^6)
# Twisted flag, g=1, d=4: (5737//2592*t1^8 - 125423//20736*t1^7*t2 + 33301//10368*t1^6*t2^2 + 194951//20736*t1^5*t2^3 - 5663//324*t1^4*t2^4 + 194951//20736*t1^3*t2^5 + 33301//10368*t1^2*t2^6 - 125423//20736*t1*t2^7 + 5737//2592*t2^8)//(t1^8 - 8*t1^7*t2 + 28*t1^6*t2^2 - 56*t1^5*t2^3 + 70*t1^4*t2^4 - 56*t1^3*t2^5 + 28*t1^2*t2^6 - 8*t1*t2^7 + t2^8)

#### Output WITHOUT experimental mode:
#
# AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}[0, (-8*t1^8 + 72*t1^7*t2 - 8*t1^7*t3 - 288*t1^6*t2^2 + 72*t1^6*t2*t3 - 8*t1^6*t3^2 + 672*t1^5*t2^3 - 288*t1^5*t2^2*t3 + 72*t1^5*t2*t3^2 - 8*t1^5*t3^3 - 1008*t1^4*t2^4 + 672*t1^4*t2^3*t3 - 288*t1^4*t2^2*t3^2 + 72*t1^4*t2*t3^3 - 8*t1^4*t3^4 + 1008*t1^3*t2^5 - 1008*t1^3*t2^4*t3 + 672*t1^3*t2^3*t3^2 - 288*t1^3*t2^2*t3^3 + 72*t1^3*t2*t3^4 - 8*t1^3*t3^5 - 672*t1^2*t2^6 + 1008*t1^2*t2^5*t3 - 1008*t1^2*t2^4*t3^2 + 672*t1^2*t2^3*t3^3 - 288*t1^2*t2^2*t3^4 + 72*t1^2*t2*t3^5 - 8*t1^2*t3^6 + 288*t1*t2^7 - 672*t1*t2^6*t3 + 1008*t1*t2^5*t3^2 - 1008*t1*t2^4*t3^3 + 672*t1*t2^3*t3^4 - 288*t1*t2^2*t3^5 + 72*t1*t2*t3^6 - 8*t1*t3^7 - 72*t2^8 + 288*t2^7*t3 - 672*t2^6*t3^2 + 1008*t2^5*t3^3 - 1008*t2^4*t3^4 + 672*t2^3*t3^5 - 288*t2^2*t3^6 + 72*t2*t3^7 - 8*t3^8)//(t1^4*t2^4 - 4*t1^4*t2^3*t3 + 6*t1^4*t2^2*t3^2 - 4*t1^4*t2*t3^3 + t1^4*t3^4 - 4*t1^3*t2^5 + 16*t1^3*t2^4*t3 - 24*t1^3*t2^3*t3^2 + 16*t1^3*t2^2*t3^3 - 4*t1^3*t2*t3^4 + 6*t1^2*t2^6 - 24*t1^2*t2^5*t3 + 36*t1^2*t2^4*t3^2 - 24*t1^2*t2^3*t3^3 + 6*t1^2*t2^2*t3^4 - 4*t1*t2^7 + 16*t1*t2^6*t3 - 24*t1*t2^5*t3^2 + 16*t1*t2^4*t3^3 - 4*t1*t2^3*t3^4 + t2^8 - 4*t2^7*t3 + 6*t2^6*t3^2 - 4*t2^5*t3^3 + t2^4*t3^4), 0, (-8//3*t1^12 + 100//3*t1^11*t2 - 4//3*t1^11*t3 - 192*t1^10*t2^2 + 52//3*t1^10*t2*t3 - 4//3*t1^10*t3^2 + 2024//3*t1^9*t2^3 - 104*t1^9*t2^2*t3 + 52//3*t1^9*t2*t3^2 - 4//3*t1^9*t3^3 - 4840//3*t1^8*t2^4 + 1144//3*t1^8*t2^3*t3 - 104*t1^8*t2^2*t3^2 + 52//3*t1^8*t2*t3^3 - 4//3*t1^8*t3^4 + 2772*t1^7*t2^5 - 2860//3*t1^7*t2^4*t3 + 1144//3*t1^7*t2^3*t3^2 - 104*t1^7*t2^2*t3^3 + 52//3*t1^7*t2*t3^4 - 4//3*t1^7*t3^5 - 3520*t1^6*t2^6 + 1716*t1^6*t2^5*t3 - 2860//3*t1^6*t2^4*t3^2 + 1144//3*t1^6*t2^3*t3^3 - 104*t1^6*t2^2*t3^4 + 52//3*t1^6*t2*t3^5 - 4//3*t1^6*t3^6 + 3344*t1^5*t2^7 - 2288*t1^5*t2^6*t3 + 1716*t1^5*t2^5*t3^2 - 2860//3*t1^5*t2^4*t3^3 + 1144//3*t1^5*t2^3*t3^4 - 104*t1^5*t2^2*t3^5 + 52//3*t1^5*t2*t3^6 - 4//3*t1^5*t3^7 - 2376*t1^4*t2^8 + 2288*t1^4*t2^7*t3 - 2288*t1^4*t2^6*t3^2 + 1716*t1^4*t2^5*t3^3 - 2860//3*t1^4*t2^4*t3^4 + 1144//3*t1^4*t2^3*t3^5 - 104*t1^4*t2^2*t3^6 + 52//3*t1^4*t2*t3^7 - 4//3*t1^4*t3^8 + 3740//3*t1^3*t2^9 - 1716*t1^3*t2^8*t3 + 2288*t1^3*t2^7*t3^2 - 2288*t1^3*t2^6*t3^3 + 1716*t1^3*t2^5*t3^4 - 2860//3*t1^3*t2^4*t3^5 + 1144//3*t1^3*t2^3*t3^6 - 104*t1^3*t2^2*t3^7 + 52//3*t1^3*t2*t3^8 - 4//3*t1^3*t3^9 - 1408//3*t1^2*t2^10 + 2860//3*t1^2*t2^9*t3 - 1716*t1^2*t2^8*t3^2 + 2288*t1^2*t2^7*t3^3 - 2288*t1^2*t2^6*t3^4 + 1716*t1^2*t2^5*t3^5 - 2860//3*t1^2*t2^4*t3^6 + 1144//3*t1^2*t2^3*t3^7 - 104*t1^2*t2^2*t3^8 + 52//3*t1^2*t2*t3^9 - 4//3*t1^2*t3^10 + 120*t1*t2^11 - 1144//3*t1*t2^10*t3 + 2860//3*t1*t2^9*t3^2 - 1716*t1*t2^8*t3^3 + 2288*t1*t2^7*t3^4 - 2288*t1*t2^6*t3^5 + 1716*t1*t2^5*t3^6 - 2860//3*t1*t2^4*t3^7 + 1144//3*t1*t2^3*t3^8 - 104*t1*t2^2*t3^9 + 52//3*t1*t2*t3^10 - 4//3*t1*t3^11 - 20*t2^12 + 120*t2^11*t3 - 1408//3*t2^10*t3^2 + 3740//3*t2^9*t3^3 - 2376*t2^8*t3^4 + 3344*t2^7*t3^5 - 3520*t2^6*t3^6 + 2772*t2^5*t3^7 - 4840//3*t2^4*t3^8 + 2024//3*t2^3*t3^9 - 192*t2^2*t3^10 + 100//3*t2*t3^11 - 8//3*t3^12)//(t1^6*t2^6 - 6*t1^6*t2^5*t3 + 15*t1^6*t2^4*t3^2 - 20*t1^6*t2^3*t3^3 + 15*t1^6*t2^2*t3^4 - 6*t1^6*t2*t3^5 + t1^6*t3^6 - 6*t1^5*t2^7 + 36*t1^5*t2^6*t3 - 90*t1^5*t2^5*t3^2 + 120*t1^5*t2^4*t3^3 - 90*t1^5*t2^3*t3^4 + 36*t1^5*t2^2*t3^5 - 6*t1^5*t2*t3^6 + 15*t1^4*t2^8 - 90*t1^4*t2^7*t3 + 225*t1^4*t2^6*t3^2 - 300*t1^4*t2^5*t3^3 + 225*t1^4*t2^4*t3^4 - 90*t1^4*t2^3*t3^5 + 15*t1^4*t2^2*t3^6 - 20*t1^3*t2^9 + 120*t1^3*t2^8*t3 - 300*t1^3*t2^7*t3^2 + 400*t1^3*t2^6*t3^3 - 300*t1^3*t2^5*t3^4 + 120*t1^3*t2^4*t3^5 - 20*t1^3*t2^3*t3^6 + 15*t1^2*t2^10 - 90*t1^2*t2^9*t3 + 225*t1^2*t2^8*t3^2 - 300*t1^2*t2^7*t3^3 + 225*t1^2*t2^6*t3^4 - 90*t1^2*t2^5*t3^5 + 15*t1^2*t2^4*t3^6 - 6*t1*t2^11 + 36*t1*t2^10*t3 - 90*t1*t2^9*t3^2 + 120*t1*t2^8*t3^3 - 90*t1*t2^7*t3^4 + 36*t1*t2^6*t3^5 - 6*t1*t2^5*t3^6 + t2^12 - 6*t2^11*t3 + 15*t2^10*t3^2 - 20*t2^9*t3^3 + 15*t2^8*t3^4 - 6*t2^7*t3^5 + t2^6*t3^6)]
# Twisted flag, g=0, d=1: 1
# Twisted flag, g=0, d=2: -7//8
# Twisted flag, g=0, d=3: 55//27
# Twisted flag, g=0, d=4: -455//64
# Twisted flag, g=1, d=1: 1//12
# Twisted flag, g=1, d=2: (-1//24*t1^4 + 5//12*t1^3*t2 - 7//8*t1^2*t2^2 + 5//12*t1*t2^3 - 1//24*t2^4)//(t1^4 - 4*t1^3*t2 + 6*t1^2*t2^2 - 4*t1*t2^3 + t2^4)
# Twisted flag, g=1, d=3: (-29//36*t1^6 + 5//2*t1^5*t2 - 29//36*t1^4*t2^2 - 13//6*t1^3*t2^3 - 29//36*t1^2*t2^4 + 5//2*t1*t2^5 - 29//36*t2^6)//(t1^6 - 6*t1^5*t2 + 15*t1^4*t2^2 - 20*t1^3*t2^3 + 15*t1^2*t2^4 - 6*t1*t2^5 + t2^6)
# Twisted flag, g=1, d=4: (499//48*t1^8 - 191//3*t1^7*t2 + 7313//48*t1^6*t2^2 - 583//3*t1^5*t2^3 + 36343//192*t1^4*t2^4 - 583//3*t1^3*t2^5 + 7313//48*t1^2*t2^6 - 191//3*t1*t2^7 + 499//48*t2^8)//(t1^8 - 8*t1^7*t2 + 28*t1^6*t2^2 - 56*t1^5*t2^3 + 70*t1^4*t2^4 - 56*t1^3*t2^5 + 28*t1^2*t2^6 - 8*t1*t2^7 + t2^8)
# Twisted flag, g=2, d=1: 1//240
# Twisted flag, g=2, d=2: (1//240*t1^6 + 1//60*t1^5*t2 - 1//6*t1^4*t2^2 + 5//16*t1^3*t2^3 - 1//6*t1^2*t2^4 + 1//60*t1*t2^5 + 1//240*t2^6)//(t1^6 - 6*t1^5*t2 + 15*t1^4*t2^2 - 20*t1^3*t2^3 + 15*t1^2*t2^4 - 6*t1*t2^5 + t2^6)
# Twisted flag, g=2, d=3: (1//48*t1^8 + 1//4*t1^7*t2 - 47//72*t1^6*t2^2 - 125//72*t1^5*t2^3 + 19//4*t1^4*t2^4 - 125//72*t1^3*t2^5 - 47//72*t1^2*t2^6 + 1//4*t1*t2^7 + 1//48*t2^8)//(t1^8 - 8*t1^7*t2 + 28*t1^6*t2^2 - 56*t1^5*t2^3 + 70*t1^4*t2^4 - 56*t1^3*t2^5 + 28*t1^2*t2^6 - 8*t1*t2^7 + t2^8)
# Twisted flag, g=2, d=4: (-289//48*t1^10 + 905//24*t1^9*t2 - 3541//36*t1^8*t2^2 + 2931//16*t1^7*t2^3 - 21151//64*t1^6*t2^4 + 249883//576*t1^5*t2^5 - 21151//64*t1^4*t2^6 + 2931//16*t1^3*t2^7 - 3541//36*t1^2*t2^8 + 905//24*t1*t2^9 - 289//48*t2^10)//(t1^10 - 10*t1^9*t2 + 45*t1^8*t2^2 - 120*t1^7*t2^3 + 210*t1^6*t2^4 - 252*t1^5*t2^5 + 210*t1^4*t2^6 - 120*t1^3*t2^7 + 45*t1^2*t2^8 - 10*t1*t2^9 + t2^10)
# Twisted flag, g=3, d=1: 1//6048
# Twisted flag, g=3, d=2: (1//864*t1^8 - 19//4320*t1^7*t2 - 77//8640*t1^6*t2^2 + 559//8640*t1^5*t2^3 - 1843//17280*t1^4*t2^4 + 559//8640*t1^3*t2^5 - 77//8640*t1^2*t2^6 - 19//4320*t1*t2^7 + 1//864*t2^8)//(t1^8 - 8*t1^7*t2 + 28*t1^6*t2^2 - 56*t1^5*t2^3 + 70*t1^4*t2^4 - 56*t1^3*t2^5 + 28*t1^2*t2^6 - 8*t1*t2^7 + t2^8)
# Twisted flag, g=3, d=3: (29//6048*t1^10 - 29//3780*t1^9*t2 - 2803//15120*t1^8*t2^2 + 11629//30240*t1^7*t2^3 + 2141//2160*t1^6*t2^4 - 11389//4320*t1^5*t2^5 + 2141//2160*t1^4*t2^6 + 11629//30240*t1^3*t2^7 - 2803//15120*t1^2*t2^8 - 29//3780*t1*t2^9 + 29//6048*t2^10)//(t1^10 - 10*t1^9*t2 + 45*t1^8*t2^2 - 120*t1^7*t2^3 + 210*t1^6*t2^4 - 252*t1^5*t2^5 + 210*t1^4*t2^6 - 120*t1^3*t2^7 + 45*t1^2*t2^8 - 10*t1*t2^9 + t2^10)


######################################################
#### END OF TESTS
#### Below this are just some previous tools used for debugging by Daniel.
######################################################

# for e in edges(G.g)
#   for d in 1:2
#     h = GKMtools._h(e, d, get_connection(G), G.equivariantCohomology, gens(G.equivariantCohomology.coeffRing), G.equivariantCohomology.edgeWeightClasses)
#     println("e = $e, d=$d => h = $(factor(numerator(h))) // $(factor(denominator(h)))")
#   end
# end

# e = Edge(2, 1), d=1 => h = -1 // 1 * (t1 - t3) * (t2 - t3) * (t1 - t2)^2
# e = Edge(2, 1), d=2 => h = 8 // 1 * (t1 - t3) * (t2 - t3) * (t1 - t2)^4 * (t1 + t2 - 2*t3)
# e = Edge(3, 1), d=1 => h = 1 // 1 * (t1 - t3)^2 * (t2 - t3) * (t1 - t2)
# e = Edge(3, 1), d=2 => h = -8 // 1 * (t1 - t3)^4 * (t2 - t3) * (t1 - t2) * (t1 - 2*t2 + t3)
# e = Edge(3, 2), d=1 => h = -1 // 1 * (t1 - t3) * (t1 - t2) * (t2 - t3)^2
# e = Edge(3, 2), d=2 => h = -4 // (1//2) * (t1 - t3) * (t1 - t2) * (t2 - t3)^4 * (2*t1 - t2 - t3)

# for Ev in 1:2
#   for Sv in 0:6
#     println("Vertex poly for Ev = $Ev, Sv = $Sv: $(factor(GKMtools.vertex_polynomial(2, Ev, Sv, 1, H)))")
#   end
# end


function test_sums(l::Vector, possible_factors::Vector)
  n = length(l)
  nf = length(possible_factors)
  for m in with_replacement_combinations(1:nf, n)
    for p in multiset_permutations(m, n)
      s = sum(l .* [possible_factors[p[i]] for i in 1:n])
      if is_zero(s)
        println("Got zero via $p")
      end
    end
  end
end

# R, (t1, t2, t3) = polynomial_ring(QQ, [:t1, :t2, :t3])

# lll = [(-32//3) * (t2 - t3)^6 // (1 * (t1 - t3) * (t1 - t2)^4 * (t1 + t2 - 2*t3)),
# (-32//3) * (t2 - t3)^5 // (1 * (t1 - t2)^4 * (t1 + t2 - 2*t3)),
# (-16//3) * (t1 - t2)^5 // ((1//2) * (t2 - t3)^4 * (2*t1 - t2 - t3)),
# (-16//3) * (t1 - t2)^6 // ((1//2) * (t1 - t3) * (t2 - t3)^4 * (2*t1 - t2 - t3)),
# (-4//3) * (t2 - t3)^5 // (1 * (t1 - t3)^2 * (t1 - t2)^3),

# (4//3) * (t2 - t3)^5 * (2*t1 - t2 - t3) // (1 * (t1 - t3)^2 * (t1 - t2)^4),
# (4//3) * (t2 - t3)^4 // (1 * (t1 - t3) * (t1 - t2)^3),

# (-4//3) * (t2 - t3)^4 * (t1 - 2*t2 + t3) // (1 * (t1 - t3) * (t1 - t2)^4),

# (1//24) * (t1 - 2*t2 + t3)^5 // (1 * (t1 - t3)^2 * (t2 - t3) * (t1 - t2)^2),
# (-1//24) * (t1 - 2*t2 + t3)^5 // (1 * (t1 - t3)^2 * (t1 - t2) * (t2 - t3)^2),
# (1//24) * (t1 - 2*t2 + t3)^6 // (1 * (t1 - t3)^2 * (t2 - t3)^2 * (t1 - t2)^2),
# (-1//48) * (t2 - t3)^4 // ((1//2) * (t1 - t3) * (t1 - t2)^2 * (2*t1 - t2 - t3)),
# (-1//48) * (t2 - t3)^4 // ((1//2) * (t1 - t3)^2 * (t1 - t2) * (2*t1 - t2 - t3)),
# (1//24) * (t2 - t3)^4 // (1 * (t1 - t3)^2 * (t1 - t2)^2),
# (-1//24) * (t1 - t2)^4 // (1 * (t1 - t3)^2 * (t2 - t3) * (t1 + t2 - 2*t3)),
# (-1//24) * (t1 - t2)^4 // (1 * (t1 - t3) * (t2 - t3)^2 * (t1 + t2 - 2*t3)),
# (1//24) * (t1 - t2)^4 // (1 * (t1 - t3)^2 * (t2 - t3)^2),
# (4//3) * (t1 - t2)^4 //( 1 * (t1 - t3) * (t2 - t3)^3),

# (4//3) * (t1 - t2)^4 * (t1 - 2*t2 + t3) // (1 * (t1 - t3) * (t2 - t3)^4),

# (-4//3) * (t1 - t2)^5 // (1 * (t1 - t3)^2 * (t2 - t3)^3),

# (4//3) * (t1 - t2)^5 * (t1 + t2 - 2*t3) // (1 * (t1 - t3)^2 * (t2 - t3)^4)]

function gw_kernel(l::Vector, maxDeg::Int64)
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  l = l .* lcm(denominator.(l))
  l = [evaluate(e, [x, zero(x), -one(x)]) for e in l]
  l = [numerator(e) * (1 // coeff(denominator(e), 1)) for e in l] # should be QQMPolyRingElem now.
  M = matrix(QQ, [coeff(e, x^i) for e in l, i in 0:maxDeg])
  return kernel(M; side=:left)
end

function gw_kernel_2var(l::Vector, maxDeg::Int64)
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  l = l .* lcm(denominator.(l))
  l = [evaluate(e, [x, one(x)]) for e in l]
  l = [numerator(e) * (1 // coeff(denominator(e), 1)) for e in l] # should be QQMPolyRingElem now.
  M = matrix(QQ, [coeff(e, x^i) for e in l, i in 0:maxDeg])
  return kernel(M; side=:left)
end

function test_F(dRange::UnitRange, gRange::UnitRange, F::GKMtools.AbstractGKM_graph, gamma::GKMtools.CurveClass_type, H::Dict{GKMtools.HodgeKey, QQFieldElem})
  for g in gRange
    for d in dRange
      println("\nGenus $g, degree $d:\n")
      gw, l = gromov_witten_pos_gen(F, d*gamma, 0, g, [class_one()], H)
      K = gw_kernel_2var(l, 100)
      println(K)
    end
  end
end

function print_h(F::Any, e::Edge, dRange::UnitRange)
  RR, (z1, z2) = polynomial_ring(QQ, 2, :z)
  con = get_any_connection(F)
  R = F.equivariantCohomology
  t = [z2, z2-z1]
  edge_weight_dict = R.edgeWeightClasses
  for d in dRange
    h = GKMtools._h(e, d, con, R, t, Dict(Edge(3, 4) => z1, Edge(4, 3) => -z1); check_degrees=true)
    println("h for $d is: $(factor(numerator(h))) // $(factor(denominator(h)))")
  end
end