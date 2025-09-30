#
# Each Hodge integral is indexed by two tuples (a_1, a_2, ..., a_n) and (l_1, ..., l_g)
# where a_1 >= a_2 >= ... >= a_n >= 0 and l_i >= 0 for all i.
# The number of marks n and the genus g are the length of those tuples, respectively.
# The resulting hodge integral is
#
# \int_{\bar{M}_{g,n}} \psi_1^{a_1} ... \psi_n^{a_n} \lambda_1^{l_1} ... \lambda_g^{l_g}
#
HodgeKey = Tuple{Tuple{Vararg{Int64}}, Tuple{Vararg{Int64}}}
_HodgeFolder = "M2/"

function hodge_integral(g::Int64, n::Int64, psi::Vector{Int64}, lambda::Vector{Int64}, H::Dict{HodgeKey, QQFieldElem})::QQFieldElem
  @req length(psi) == n "psi must have length n"
  @req length(lambda) == g "lambda must have length g"
  if sum(psi) + sum((1:g) .* lambda) != 3*g - 3 + n
    return zero(QQ)
  end
  return H[((sort(psi; rev=true)...,), (lambda...,))]
end

# Desired format for storing Hodge integrals:
#
# g, n, (psi (descending)), (l1, l2, ..., lg); hodgeIntegral

function _load_Hodge_integrals(dir::String)::Dict{HodgeKey, QQFieldElem}

  H = Dict{HodgeKey, QQFieldElem}()

  for f in readdir(dir)
    filename = dir * f
    # println("Loading Hodge integrals from $filename")
    for l in eachline(filename)
      _add_to_Hodge_dict(l, H)
    end
  end
  return H
end

function _add_to_Hodge_dict(line::String, H::Dict{HodgeKey, QQFieldElem})
  S = split(line, ";")
  @req length(S) == 2 "Line has wrong format."
  ind, val = S

  # parse value of Hodge integral.
  val = split(val, "//")
  @req length(val) == 2 "Invalid format after ;"
  val = ZZ(val[1]) // ZZ(val[2])

  # parse index of Hodge integral.
  ind = split(ind, ",")
  @req length(ind) >= 4 "Invalid format before ;"
  g = parse(Int64, ind[1])
  n = parse(Int64, ind[2])
  for i in 1:length(ind)
    ind[i] = replace(ind[i], "(" => "")
    ind[i] = replace(ind[i], ")" => "")
    ind[i] = replace(ind[i], " " => "")
  end
  filter!(i -> i != "", ind)
  @req length(ind) == 2 + g + n "Wrong number of integers before ;"
  psi = parse.(Int64, ind[3:2+n])
  psi = (psi...,)
  lambda = parse.(Int64, ind[3+n:2+n+g])
  lambda = (lambda...,)

  key = (psi, lambda)
  if haskey(H, key)
    @req H[key] == val "Conflicting values found for index $key"
  end

  H[key] = val
end