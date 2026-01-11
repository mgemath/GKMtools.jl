# Entry [g, n] is the genus g vertex poly with n marked points and without psi classes,
# for the fixed value of valG.
function vertex_polynomials(gMax::Int64, nMax::Int64, valG::Int64)::Matrix{QQMPolyRingElem}
  @req gMax >= 0 "gMax must be non-negative."
  @req nMax > 0 "nMax must be positive"
  @req valG > 0 "valG must be positive"

  H = load_H(gMax, nMax)
  return vertex_polynomials(gMax, nMax, valG, H)
end

function vertex_polynomials(gMax::Int64, nMax::Int64, valG::Int64, H::Dict{HodgeKey, QQFieldElem})::Matrix{QQMPolyRingElem}
  @req gMax >= 0 "gMax must be non-negative."
  @req nMax > 0 "nMax must be positive"
  @req valG > 0 "valG must be positive"

  VPs = Matrix{QQMPolyRingElem}(undef, gMax, nMax)
  R, w, u = polynomial_ring(QQ, ["w$i" for i in 1:valG], ["u$i" for i in 1:nMax])
  for g in 1:gMax
    for n in 1:nMax
      VPs[g, n] = _vertex_polynomial(valG, n, zeros(Int64, 0), zeros(Int64, g), g, H, w, u; prefactor=true)
    end
  end
  return VPs
end

# This is called by Euler_inv_pos_gen(...).
# Use this version of the function to evaluate a vertex polynomial using the output of
# vertex_polynomials(...).
function evaluate_vertex_polynomial(u::Vector{T}, w::Vector{T}, nMarks::Int64, g::Int64, VPs::Matrix{QQMPolyRingElem}) where T<:RingElem
  # valG = length(w)
  Ev = length(u)
  res = zero(w[1])
  if iszero(g)
    res = prod(u) * (Ev-3+nMarks >= 0 ? sum(u)^(Ev-3+nMarks) : 1//sum(u)^(-(Ev-3+nMarks)))
  else
    vp = VPs[g, Ev]
    _, nMax = size(VPs) # gMax, nMax = size(vp)
    res = evaluate(vp, vcat(w, u, repeat([zero(w[1])], nMax - Ev)))
    res *= sum(u)^nMarks # String equation
  end
  # println("nMarks=$nMarks, value= ", res)
  return res
end

# This is called by Euler_inv_pos_gen_2(...).
function evaluate_vertex_polynomial_with_psis(u::Vector{T}, w::Vector{T}, psi_exp::Vector{Int64}, g::Int64, H::Dict{HodgeKey, QQFieldElem}) where T<:RingElem
  Ev = length(u)
  Sv = length(psi_exp)
  n = Ev+Sv
  res = zero(w[1])
  if iszero(g)
    sum_psi = sum(psi_exp)
    e = n - 3 - sum_psi
    combinatorial_factor = prod(e+1:n-3) // (Sv == 0 ? 1 : prod(a -> factorial(a), psi_exp))
    res = combinatorial_factor * prod(u) * (e >= 0 ? sum(u)^(e) : 1//sum(u)^(-e))
  else
    # vp = VPs[g, Ev]
    # _, nMax = size(VPs) # gMax, nMax = size(vp)
    # res = evaluate(vp, vcat(w, u, repeat([zero(w[1])], nMax - Ev)))
    valG = length(w)
    res = _vertex_polynomial(valG, Ev, psi_exp, zeros(Int64, g), g, H, w, u)
  end
  # println("nMarks=$nMarks, value= ", res)
  return res
end

# To evaluate a single vertex polynomial directly form the Hodge numbers, use this function.
function evaluate_vertex_polynomial(u::Vector{T}, w::Vector{T}, nMarks::Int64, g::Int64, H::Dict{HodgeKey, QQFieldElem}) where T<:RingElem
  vp = vertex_polynomial(length(w), length(u), nMarks, g, H; prefactor=true)
  # println("nMarks=$nMarks, vp= $vp, value= ", evaluate(vp, vcat(w, u)))
  return evaluate(vp, vcat(w, u))
end

function vertex_polynomial(valG::Int64, Ev::Int64, Sv::Int64, gv::Int64, H::Dict{HodgeKey, QQFieldElem}; prefactor::Bool=true)
  @req Sv >= 0 "Sv must be non-negative"
  return vertex_polynomial(valG, Ev, zeros(Int64, Sv), zeros(Int64, gv), gv, H; prefactor=prefactor)
end

# To get something homogeneous, the result is divided by w_\epsilon^{g} for each \epsilon, and
# the substitution w -> 1/w is applied.
function vertex_polynomial(valG::Int64, Ev::Int64, markPsis::Vector{Int64}, lambda::Vector{Int64}, gv::Int64, H::Dict{HodgeKey, QQFieldElem}; prefactor::Bool=true)
  @req all(i -> i >= 0, markPsis) "markPsi entries must be non-negative"
  @req all(i -> i >= 0, lambda) "all lambda entries must be non-negative"
  @req gv >= 0 "gv must be non-negative"
  @req Ev >= 1 "Ev must be positive"
  @req valG >= 1 "valG must be positive"

  _, w, u = polynomial_ring(QQ, ["w$i" for i in 1:valG], ["u$i" for i in 1:Ev])
  return _vertex_polynomial(valG, Ev, markPsis, lambda, gv, H, w, u; prefactor=prefactor)
end

# To get something homogeneous, the result is divided by w_\epsilon^{g} for each \epsilon, and
# the substitution w -> 1/w is applied.
# supports lambda insertions as well now.
function _vertex_polynomial(valG::Int64, Ev::Int64, markPsis::Vector{Int64}, lambda::Vector{Int64}, gv::Int64, H::Dict{HodgeKey, QQFieldElem}, w::Vector{T}, u::Vector{T}; prefactor::Bool=true) where T<:RingElem

  Sv = length(markPsis)
  @req all(i -> i >= 0, markPsis) "markPsi entries must be non-negative"
  @req gv >= 0 "gv must be non-negative"
  @req Ev >= 1 "Ev must be positive"
  @req valG >= 1 "valG must be positive"
  @req length(lambda) == gv "length(lambda) must be gv"

  R = parent(w[1])

  res = zero(R)

  n = Ev + Sv
  g = gv
  dimM = 3*g - 3 + n
  # psi_buf = Vector{Int}(undef, n)   # n = Ev + length(markPsis)
  # lambda_buf = similar(lambda)   # length g

  # In genus zero, Liu--Sheshmani formaly allow dimM < 0. We implement this exception here (when there are no psi classes).
  if g == 0 && all(x -> iszero(x), markPsis) && all(x -> iszero(x), lambda) && dimM < 0
    #println("Using exception for Ev=$Ev, Sv=$Sv, gv=$gv")
    return (prefactor ? prod(u[1:Ev]) : one(R)) // (sum(u[1:Ev])^(-dimM))
  end

  @req dimM >= 0 "dimM must be non-negative. Got Ev=$Ev, Sv=$Sv, gv=$gv"
  totalPsiMarksAndLambda = sum(markPsis) + sum(collect(1:gv) .* lambda)
  #@req dimM - totalPsiMarks >= 0 "Too many psi marks."
  dimM - totalPsiMarksAndLambda < 0 && return zero(R)

  # C[1:valG] are exponents of w1, w2, ...
  # C[valG+1:valG+Ev] .+ 1 are exponents of u1, u2, ...
  for C in weak_compositions(dimM - totalPsiMarksAndLambda, valG + Ev)

    ##############################
    ### Optimized version with no allocations

    # any(i -> C[i] > g, 1:valG) && continue

    # @inbounds for j in 1:g
    #     lambda_buf[j] = lambda[j]
    # end

    # @inbounds for k in 1:valG
    #   lambda_buf[C[k]] += 1
    # end
    
    # m = prod( i -> w[i]^(C[i]),1:valG) * prod(i-> u[i]^(C[valG+i] + (prefactor ? 1 : 0)), 1:Ev )
    # res += hodge_integral(g, n, C, valG, Ev, markPsis, lambda_buf, H, psi_buf) * m * QQ(-1)^(sum(i -> C[i], 1:valG))
##############################

    any(i -> i > g, C[1:valG]) && continue
    
    l = C[1:valG]
    lambda_for_H = [count(i -> i==j, l) for j in 1:g] .+ lambda
    psi = vcat(C[valG+1:valG+Ev], markPsis)

    #m = prod( w.^(g .- C[1:valG]) ) * prod( u.^(C[valG+1:valG+Ev] .+ 1) )
    m = prod( w[1:valG].^(C[1:valG]) ) * prod( u[1:Ev].^(C[valG+1:valG+Ev] .+ (prefactor ? 1 : 0)) )

    res += hodge_integral(g, n, psi, lambda_for_H, H) * m * QQ(-1)^(sum(l))
  end

  return res
end

#########################
#####
#####   Below are some tools used during development to validate and experiment with Hodge integral data.
#####   Nothing below is used in the package.
#####
#########################


# The outcome of this function is:
# 1. For markPsis = [0, ..., 0] we always get u1 + u2 + ...
# 2. For g > 0 and some markPsis non-zero, this does not hold.
# 3. For g=0 it is still u1+u2+... times some factor (as we already know.)
function _test_step(valG::Int64, Ev::Int64, markPsis::Vector{Int64}, gv::Int64, H::Dict{HodgeKey, QQFieldElem})
  return vertex_polynomial(valG, Ev, vcat(markPsis, [0]), gv, H) // vertex_polynomial(valG, Ev, markPsis, gv, H) 
end


#
# The following applies only to the special case when the flags of the vertex v on the tree \Gamma 
# are mapped bijectively to the flags of \sigma_v on the GKM graph G, and all edge multiplicities are 1.
#
function _vertex_polynomial_old(gv::Int64, Ev::Int64, Sv::Int64, H::Dict{HodgeKey, QQFieldElem})
  @req gv >= 0 "gv must be non-negative"
  @req Ev >= 1 "Ev must be positive"
  @req Sv >= 0 "Sv must be non-negative"

  R, w = polynomial_ring(QQ, Ev, :w)
  res = zero(R)

  n = Ev + Sv
  g = gv
  dimM = 3*g - 3 + n
  @req dimM >= 0 "dimM must be non-negative."

  for C in weak_compositions(dimM, Ev)
    tmp = QQ(0)
    for l in Iterators.product([0:g for i in 1:Ev]...)
      psi = vcat(C .- l, zeros(Int64, Sv))
      any(i -> i < 0, psi) && continue
      lambda = [count(i -> i==j, l) for j in 1:g]
      tmp += hodge_integral(g, n, psi, lambda, H) * QQ(-1)^sum(l)
    end
    m = one(R)
    for i in 1:Ev
      m *= w[i]^(C[i])
    end
    res += tmp * m
  end
  #if g > 0
  #  res = res // prod(w[1:Ev])^(g-1)
  #  if denominator(res) == 1
  #    res = numerator(res) #
  #  end
  #else 
  #  res = res * prod(w[1:Ev])
  #end
  return res
end

# Note that the resulting polynomial is not always divisible by (w_1...w_Ev), e.g.:
#
# julia> GKMtools.vertex_polynomial(2, 3, 0, H)
# (1//1152*w1^4*w2^2 + 1//576*w1^4*w2*w3 + 1//1152*w1^4*w3^2 + 1//576*w1^3*w2^3 + 1//120*w1^3*w2^2*w3 + 1//120*w1^3*w2*w3^2 + 1//576*w1^3*w3^3 + 1//1152*w1^2*w2^4 + 1//120*w1^2*w2^3*w3 + 17//960*w1^2*w2^2*w3^2 + 1//120*w1^2*w2*w3^3 + 1//1152*w1^2*w3^4 + 1//576*w1*w2^4*w3 + 1//120*w1*w2^3*w3^2 + 1//120*w1*w2^2*w3^3 + 1//576*w1*w2*w3^4 + 1//1152*w2^4*w3^2 + 1//576*w2^3*w3^3 + 1//1152*w2^2*w3^4)//(w1*w2*w3)

function _test_hodge_divisibility(H::Dict{HodgeKey, QQFieldElem})
  A277001 = [1, 24, 5760, 2903040, 1393459200, 367873228800]
  for h in H
    k = h[1]
    psi = k[1]
    lambda = k[2]
    g = length(lambda)
    n = length(psi)
    val = h[2]
    if denominator(val * A277001[g+1]) != 1
      println(h)
      println(" ... denominator not divisible by $(A277001[g+1])")
    end
  end
end

function _hodge_denominators_by_genus(H::Dict{HodgeKey, QQFieldElem})
  denoms = Dict{Int64, Int64}()
  for h in H
    k = h[1]
    psi = k[1]
    lambda = k[2]
    g = length(lambda)
    n = length(psi)
    val = h[2]
    if !haskey(denoms, g)
      denoms[g] = 1
    end
    newVal = lcm(denoms[g], denominator(val))
    if denoms[g] != newVal
      println("Update g=$g from $(denoms[g]) to $newVal for reason:")
      println("  $h")
    end
    denoms[g] = lcm(denoms[g], denominator(val))
  end
  return sort(denoms)
end

##########
# The following function tests the string equation in our dataset of hodge integrals.
# Output: 
##########
#
# String equation fails for:
#   Pair{Tuple{Tuple{Vararg{Int64}}, Tuple{Vararg{Int64}}}, QQFieldElem}(((0, 0, 0), ()), 1)
#  sum is: 0
# String equation fails for:
#   Pair{Tuple{Tuple{Vararg{Int64}}, Tuple{Vararg{Int64}}}, QQFieldElem}(((0,), (1,)), 1//24)
#  sum is: 0
# String equation holds for 17397 out of 17399 entries with zero psi.
#
##########
#
# Explanation:
#    The first case ((0, 0, 0), ()) is relevant and must be coded as exception. (see Witten's 1990 paper.)
#    The second case ((0,), (1,)) does not happen if we assume |Ev| >= 1, so we will not 
function _check_string_equation(H::Dict{HodgeKey, QQFieldElem})
  ctr = 0
  ctr2 = 0
  for h in H
    k = h[1]
    val = h[2]
    psi = k[1]
    lambda = k[2]
    n = length(psi)
    g = length(lambda)
    if iszero(psi[n])
      ctr2 += 1
      s = QQ(0)
      for i in 1:n-1
        iszero(psi[i]) && continue
        psi_minus_i = [psi[j] for j in 1:n-1]
        psi_minus_i[i] -= 1
        s += hodge_integral(g, n-1, psi_minus_i, Int64[lambda...], H)
      end
      if s != val
        println("String equation fails for:\n  $h")
        println("  sum is: $s")
      else 
        ctr += 1
      end
    end
  end
  println("String equation holds for $ctr out of $ctr2 entries with zero psi.")
end

function _check_string_equation_for_vertex_polys(valG, Ev, gv, H::Dict{HodgeKey, QQFieldElem}; print::Bool=true)
  p = vertex_polynomial(valG, Ev, 0, gv, H; prefactor=false)
  q = vertex_polynomial(valG, Ev+1, 0, gv, H; prefactor=false)
  g = gens(parent(p))
  q0 = evaluate(q, vcat(g, [zero(g[1])]))
  s = sum(g[valG+1:valG+Ev])
  f = (s * p) == q0

  print && println("s = $s")
  print && println("p = $(factor(p))")
  print && println("q = $(factor(q))")
  print && println("q0 = $(iszero(q0) ? 0 : factor(q0)))")
  print && println(f)

  return f
end

# julia> GKMtools._check_string_equation_for_all_vertex_polys(4, 3, 5, H)
# Polynomial string equation holds for 60 out of 60 checked values.
function _check_string_equation_for_all_vertex_polys(valG_max, Ev_max, gv_max, H::Dict{HodgeKey, QQFieldElem})
  ctr = 0
  tot = 0
  for valG in 1:valG_max
    for Ev in 1:Ev_max
      for gv in 1:gv_max
        tot += 1
        if _check_string_equation_for_vertex_polys(valG, Ev, gv, H; print=false)
          ctr += 1
          continue
        end
        println("Poly string equation fails for valG=$valG, Ev=$Ev, gv=$gv:")
        _check_string_equation_for_vertex_polys(valG, Ev, gv, H; print=true)
      end
    end
  end
  println("Polynomial string equation holds for $ctr out of $tot checked values.")
end


# This function was used to test evaluate_vertex_polynomial_with_psis().
# all tests pass for input (4, 4, 4, load_H()).
function _check_genus_zero(valG_max, Ev_max, Sv_max, H::Dict{HodgeKey, QQFieldElem})
  a_max = 10
  for valG in 1:valG_max
    for Ev in 1: Ev_max
      for Sv in 1:Sv_max
        Ev + Sv < 3 && continue
        for psi_exp in Combinatorics.with_replacement_combinations(1:a_max, Sv)
          println("valG=$valG, Ev=$Ev, Sv=$Sv, psi_exp=$psi_exp")
          vp1 = vertex_polynomial(valG, Ev, psi_exp, zeros(Int64, 0), 0, H)
          w = gens(parent(vp1))[1:valG]
          u = gens(parent(vp1))[valG+1:valG+Ev]
          vp2 = evaluate_vertex_polynomial_with_psis(u, w, psi_exp, 0, H)
          println("vp1: $vp1")
          println("vp2: $vp2")
          println("equal: $(vp1 == vp2)")
          if vp1 != vp2
            error("Not equal for valG=$valG, Ev=$Ev, Sv=$Sv, psi_exp=$psi_exp")
          end
        end
      end
    end
  end
  println("All tests pass")
end

# This function was used to test evaluate_vertex_polynomial_with_psis().
# all tests pass for input (3, 3, 3, 3, load_H()).
function _check_positive_genus(max_g, valG_max, Ev_max, Sv_max, H::Dict{HodgeKey, QQFieldElem})
  a_max = 10
  for g in 1:max_g
    for valG in 1:valG_max
      for Ev in 1: Ev_max
        for Sv in 1:Sv_max
          3*g + Ev + Sv - 3 < 0 && continue
          for psi_exp in Combinatorics.with_replacement_combinations(1:a_max, Sv)
            println("g=$g, valG=$valG, Ev=$Ev, Sv=$Sv, psi_exp=$psi_exp")
            vp1 = vertex_polynomial(valG, Ev, psi_exp, zeros(Int64, g), g, H)
            w = gens(parent(vp1))[1:valG]
            u = gens(parent(vp1))[valG+1:valG+Ev]
            vp2 = evaluate_vertex_polynomial_with_psis(u, w, psi_exp, g, H)
            println("vp1: $vp1")
            println("vp2: $vp2")
            println("equal: $(vp1 == vp2)")
            if vp1 != vp2
              error("Not equal for g=$g, valG=$valG, Ev=$Ev, Sv=$Sv, psi_exp=$psi_exp")
            end
          end
        end
      end
    end
  end
  println("All tests pass")
end