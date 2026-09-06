# Mihalcea, arXiv:math/0501213, Cor. 6.5 and Sections 7-8.
# Same-degree classical steps increase l(u)+l(v)-l(w); quantum steps
# decrease the effective multidegree. Never set equivariant parameters to zero.
function _qs_recur(Q, u, v, w, d, getcoefficient)
  u == w && ((u,v) = (v,u))
  u != w || error("diagonal coefficient requires unit reconstruction")
  value = getcoefficient(0, 0, 0, d) # typed additive zero
  for e in Q.outgoing[u]
    dd = Tuple(d[i]-e.degree[i] for i in eachindex(d))
    all(>=(0), dd) || continue
    value += sum(e.coefficients)*getcoefficient(e.vertex, v, w, dd)
  end
  for e in Q.incoming[w]
    dd = Tuple(d[i]-e.degree[i] for i in eachindex(d))
    all(>=(0), dd) || continue
    value -= sum(e.coefficients)*getcoefficient(u, v, e.vertex, dd)
  end
  return value / (Q.diagonal[w]-Q.diagonal[u])
end

_qs_below(Q,u,w) = u == w || Q.representatives[u] < Q.representatives[w]

function _qs_coefficient(Q, u, v, w, d)
  u == 0 && return QQ(0)
  r = Q.codimensions[u]+Q.codimensions[v]-Q.codimensions[w]-_qs_degree_value(Q,d)
  r < 0 && return QQ(0)
  if all(iszero,d)
    (_qs_below(Q,u,w) && _qs_below(Q,v,w)) || return QQ(0)
  else
    # Mihalcea's main vanishing lemma, applied to either input.
    max(Q.codimensions[u],Q.codimensions[v])+1 > Q.codimensions[w]+_qs_degree_value(Q,d) && return QQ(0)
  end
  if u == Q.unit || v == Q.unit
    return QQ(all(iszero,d) && (u == Q.unit ? v : u) == w)
  end
  u > v && ((u,v) = (v,u))
  key = (u,v,w,d)
  haskey(Q.cache,key) && return Q.cache[key]
  value = if u == v == w
    _qs_diagonal(Q,w,d)
  else
    _qs_recur(Q,u,v,w,d,(a,b,c,e) -> _qs_coefficient(Q,a,b,c,e))
  end
  Q.cache[key] = value
  return value
end

# Represent the unit recurrence as A*x+B, x=C[w,w,w,d].
# Polynomial pairs are local to this diagonal solve. At the same degree no
# vanishing-by-polynomial-degree shortcut is allowed: the unit equation itself
# has negative degree for d>0, but determines x by cancellation.
function _qs_diagonal(Q, target, d)
  P, x = polynomial_ring(QQ, "_diagonal")
  memo = Dict{Tuple{Int,Int,Int},typeof(x)}()
  function affine(u,v,w,dd)
    u == 0 && return zero(P)
    dd != d && return P(_qs_coefficient(Q,u,v,w,dd))
    if all(iszero,d)
      (_qs_below(Q,u,w) && _qs_below(Q,v,w)) || return zero(P)
    end
    u == v == w == target && return x
    u > v && ((u,v) = (v,u))
    key = (u,v,w)
    haskey(memo,key) && return memo[key]
    result = if u == v == w
      P(_qs_coefficient(Q,u,v,w,dd))
    elseif u == Q.unit && v != w || v == Q.unit && u != w
      zero(P)
    else
      _qs_recur(Q,u,v,w,dd,affine)
    end
    memo[key] = result
    return result
  end
  expression = _qs_recur(Q,Q.unit,target,target,d,affine)
  a = coeff(expression,1)
  iszero(a) && error("singular unit reconstruction")
  rhs = all(iszero,d) ? QQ(1) : QQ(0)
  return (rhs-coeff(expression,0))/a
end

"""
    quantum_schubert_coefficient(Q, u, v, w, degree=0)

Return the ordinary coefficient of `q^degree sigma_w` in
`sigma_u ⋆ sigma_v` as an exact rational integer. Indices are context vertex
indices or representative strings. Degrees are nonnegative tuples ordered by
`Q.omitted_roots`; an integer is accepted for Picard rank one (zero for any rank).
"""
function quantum_schubert_coefficient(Q::QuantumSchubertContext,u,v,w,degree=0)
  u,v,w = _qs_vertex(Q,u),_qs_vertex(Q,v),_qs_vertex(Q,w)
  d = _qs_degree(Q,degree)
  Q.codimensions[u]+Q.codimensions[v] == Q.codimensions[w]+_qs_degree_value(Q,d) || return QQ(0)
  if Q.codimensions[u] == 1 || Q.codimensions[v] == 1
    a,b = Q.codimensions[u] == 1 ? (u,v) : (v,u)
    i = Int(only(word(Q.representatives[a])))
    return get(quantum_chevalley_product(Q,i,b),(w,d),QQ(0))
  end
  c = _qs_coefficient(Q,u,v,w,d)
  denominator(c) == 1 && c >= 0 || error("reconstruction produced a nonintegral or negative ordinary coefficient")
  return c
end

"""
    quantum_schubert_product(Q, u, v; max_degree=nothing)

Compute an ordinary quantum basis product as a sparse dictionary
`(output_vertex, degree_tuple) => QQ`. By default all possible degrees are
included. `max_degree` is an optional componentwise Novikov bound.
"""
function quantum_schubert_product(Q::QuantumSchubertContext,u,v;max_degree=nothing)
  u,v = _qs_vertex(Q,u),_qs_vertex(Q,v)
  cap = isnothing(max_degree) ? nothing : _qs_degree(Q,max_degree)
  if Q.codimensions[u] == 1 || Q.codimensions[v] == 1
    a,b = Q.codimensions[u] == 1 ? (u,v) : (v,u)
    i = Int(only(word(Q.representatives[a])))
    result = quantum_chevalley_product(Q,i,b)
    isnothing(cap) || filter!(p -> all(first(p)[2][j] <= cap[j] for j in eachindex(cap)), result)
    return result
  end
  result = Dict{Tuple{Int,_QSDegree},QQFieldElem}()
  total = Q.codimensions[u]+Q.codimensions[v]
  for d in _qs_degrees(Q,total)
    isnothing(cap) || all(d[i] <= cap[i] for i in eachindex(d)) || continue
    outdegree = total-_qs_degree_value(Q,d)
    for w in eachindex(Q.representatives)
      Q.codimensions[w] == outdegree || continue
      c = quantum_schubert_coefficient(Q,u,v,w,d)
      iszero(c) || (result[(w,d)] = c)
    end
  end
  return result
end

"""
    quantum_schubert_table(Q; max_degree=nothing)

Compute every ordinary basis product. The dictionary stores pairs `(u,v)`
with `u <= v`, each mapped to the sparse output of
[`quantum_schubert_product`](@ref). This can be large; individual products
share the same lazy reconstruction cache.
"""
function quantum_schubert_table(Q::QuantumSchubertContext;max_degree=nothing)
  return Dict((u,v) => quantum_schubert_product(Q,u,v;max_degree)
    for u in eachindex(Q.representatives) for v in u:length(Q.representatives))
end
