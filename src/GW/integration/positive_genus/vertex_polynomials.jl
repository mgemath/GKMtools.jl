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
  totalPsiMarksAndLambda = sum(markPsis) + sum(i * lambda[i] for i in eachindex(lambda))
  #@req dimM - totalPsiMarks >= 0 "Too many psi marks."
  dimM - totalPsiMarksAndLambda < 0 && return zero(R)

  lambda_buffer = similar(lambda)
  psi_buffer = Vector{Int}(undef, n)

  # C[1:valG] are exponents of w1, w2, ...
  # C[valG+1:valG+Ev] are exponents of u1, u2, ...
  for C in weak_compositions(dimM - totalPsiMarksAndLambda, valG + Ev)

    copyto!(lambda_buffer, lambda)
    valid = true
    sign_exponent = 0
     for i in 1:valG
      exponent = C[i]
      if exponent > g
        valid = false
        break
      end
      exponent > 0 && (lambda_buffer[exponent] += 1)
      sign_exponent += exponent
    end
    valid || continue

     for i in 1:Ev
      exponent = C[valG + i]
      psi_buffer[i] = exponent
    end
     for i in eachindex(markPsis)
      psi_buffer[Ev + i] = markPsis[i]
    end
    sort!(psi_buffer; rev=true)
    hodge = get(H, ((psi_buffer...,), (lambda_buffer...,)), zero(QQ))
    iszero(hodge) && continue

    # Construct the ring monomial only for a nonzero Hodge coefficient.
    m = one(R)
    for i in 1:valG
      exponent = C[i]
      exponent > 0 && (m *= w[i]^exponent)
    end
    for i in 1:Ev
      exponent = C[valG + i]
      u_exponent = exponent + (prefactor ? 1 : 0)
      u_exponent > 0 && (m *= u[i]^u_exponent)
    end
    res += (isodd(sign_exponent) ? -hodge : hodge) * m
  end

  return res
end
