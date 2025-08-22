# Return BPS state using Gromov--Witten invariants and Gopakumar--Vafa formula.
function BPS_states(k::Int64, dMax::Int64; printStuff::Bool = true)
  Nd = [(-1)^( (k+1)*d + 1 ) // (k^2 * d^3) * binomial(k^2 * d, d) for d in 1:dMax]
  nd = Vector(undef, dMax)
  for d in 1:dMax
    nd[d] = Nd[d]
    for e in 1:(d-1)
      !(divides(d, e)[1]) && continue
      # e divides d and e < d.
      nd[d] -= 1//(d // e)^3 * nd[e]
    end
    if denominator(nd[d] == 1)
      nd[d] = numerator(nd[d])
    end
  end
  printStuff && println("k = $k")
  printStuff && println("Nd = $Nd")
  printStuff && println("nd = $nd")
  return nd
end

# Return BPS state using Mobius inversion.
function BPS_state_mobius(k::Int64, d::Int64)
  r = 0
  for e in 1:d
    !(divides(d, e)[1]) && continue
      # e divides d and e < d.
    m = mobius(e)
    de = Int64(d//e)
    m == 0 && continue
    r += m * (-1)^( (k+1)*de + 1 ) * binomial(k^2 * de, de)
  end
  r //= d^3 * k^2
  if denominator(r) == 1
    return numerator(r)
  end
  return r
end

# Return the Mobius function evaluated at n.
function mobius(n::Int64)
  m = 1
  for f in factor(n)
    f[2] > 1 && return 0
    m *= -1
  end
  return m
end

# Print the Table 5.1.
function print_BPS_table(kMax::Int64, dMax::Int64; latex::Bool=false)
  for k in 1:kMax
    nd = BPS_states(k, dMax; printStuff=false)
    for d in 1:dMax
      @req nd[d] == BPS_state_mobius(k, d) "BPS states disagree for d=$d, k=$k"
    end
    if latex
      print("\$$k\$ & ")
      for d in 1:dMax
        print("\$$(nd[d])\$")
        d != dMax && print(" & ")
      end
      println(" \\\\\\hline")
    else
      print("k = $k:\t")
      for d in 1:dMax
        print("$(nd[d])\t")
      end
      println()
    end
  end
end