# This is called by Euler_inv_pos_gen(...).
function evaluate_vertex_polynomial(u::Vector{T}, w::Vector{T}, nMarks::Int64, g::Int64, H::Dict{HodgeKey, QQFieldElem}) where T<:RingElem
  vp = vertex_polynomial(length(w), length(u), nMarks, g, H; prefactor=true)
  return evaluate(vp, vcat(w, u))
end

function vertex_polynomial(valG::Int64, Ev::Int64, Sv::Int64, gv::Int64, H::Dict{HodgeKey, QQFieldElem}; prefactor::Bool=true)
  @req Sv >= 0 "Sv must be non-negative"
  return vertex_polynomial(valG, Ev, zeros(Int64, Sv), gv, H; prefactor=prefactor)
end

# To get something homogeneous, the result is divided by w_\epsilon^{g} for each \epsilon, and
# the substitution w -> 1/w is applied.
function vertex_polynomial(valG::Int64, Ev::Int64, markPsis::Vector{Int64}, gv::Int64, H::Dict{HodgeKey, QQFieldElem}; prefactor::Bool=true)

  Sv = length(markPsis)
  @req all(i -> i >= 0, markPsis) "markPsi entries must be non-negative"
  @req gv >= 0 "gv must be non-negative"
  @req Ev >= 1 "Ev must be positive"
  @req valG >= 1 "valG must be positive"

  R, w, u = polynomial_ring(QQ, ["w$i" for i in 1:valG], ["u$i" for i in 1:Ev])
  res = zero(R)

  n = Ev + Sv
  g = gv
  dimM = 3*g - 3 + n

  # In genus zero, Liu--Sheshmani formaly allow dimM < 0. We implement this exception here (when there are no psi classes).
  if g == 0 && all(x -> iszero(x), markPsis) && dimM < 0
    #println("Using exception for Ev=$Ev, Sv=$Sv, gv=$gv")
    return (prefactor ? prod(u) : one(R)) // (sum(u)^(-dimM))
  end

  @req dimM >= 0 "dimM must be non-negative. Got Ev=$Ev, Sv=$Sv, gv=$gv"
  totalPsiMarks = sum(markPsis)
  @req dimM - totalPsiMarks >= 0 "Too many psi marks."

  # C[1:valG] are exponents of w1, w2, ...
  # C[valG+1:valG+Ev] .+ 1 are exponents of u1, u2, ...
  for C in weak_compositions(dimM - totalPsiMarks, valG + Ev)
    any(i -> i > g, C[1:valG]) && continue
    
    l = C[1:valG]
    lambda = [count(i -> i==j, l) for j in 1:g]
    psi = vcat(C[valG+1:valG+Ev], markPsis)

    #m = prod( w.^(g .- C[1:valG]) ) * prod( u.^(C[valG+1:valG+Ev] .+ 1) )
    m = prod( w.^(C[1:valG]) ) * prod( u.^(C[valG+1:valG+Ev] .+ (prefactor ? 1 : 0)) )

    #println(m)
    #println(hodge_integral(g, n, psi, lambda, H) * m * QQ(-1)^(sum(l)))

    #println("psi = $psi, lambda=$lambda")
    res += hodge_integral(g, n, psi, lambda, H) * m * QQ(-1)^(sum(l))
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
function vertex_polynomial_old(gv::Int64, Ev::Int64, Sv::Int64, H::Dict{HodgeKey, QQFieldElem})
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