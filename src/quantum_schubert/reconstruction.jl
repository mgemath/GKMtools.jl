# Mihalcea, arXiv:math/0501213, Cor. 6.5 and Sections 7-8.
# Same-degree classical steps increase l(u)+l(v)-l(w); quantum steps
# decrease the effective multidegree. Never set equivariant parameters to zero.
function _qs_recur(Q, u, v, w, d, getcoefficient)
  # Start from the higher-codimension input: fewer classical upward steps
  # remain before the vanishing bound. Canonical cache keys are independent
  # of this choice of associativity equation.
  (u == w || (v != w && Q.codimensions[v] > Q.codimensions[u])) && ((u,v) = (v,u))
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

_qs_below(Q,u,w) = isnothing(Q.lower_intervals) ? (u == w || Q.representatives[u] < Q.representatives[w]) : Q.lower_intervals[w][u]

function _qs_coefficient(Q, u, v, w, d)
  u == 0 && return zero(Q.field)
  r = Q.codimensions[u]+Q.codimensions[v]-Q.codimensions[w]-_qs_degree_value(Q,d)
  r < 0 && return zero(Q.field)
  if all(iszero,d)
    (_qs_below(Q,u,w) && _qs_below(Q,v,w)) || return zero(Q.field)
    u == w && return _qs_restriction(Q,v,w)
    v == w && return _qs_restriction(Q,u,w)
  else
    # Mihalcea's main vanishing lemma, applied to either input.
    max(Q.codimensions[u],Q.codimensions[v])+1 > Q.codimensions[w]+_qs_degree_value(Q,d) && return zero(Q.field)
  end
  if u == Q.unit || v == Q.unit
    return Q.field(all(iszero,d) && (u == Q.unit ? v : u) == w)
  end
  u > v && ((u,v) = (v,u))
  key = (u,v,w,d)
  haskey(Q.cache,key) && return Q.cache[key]
  value = if u == v == w
    _qs_diagonal(Q,w,d)
  else
    _qs_recur(Q,u,v,w,d,(a,b,c,e) -> _qs_coefficient(Q,a,b,c,e))
  end
  return _qs_remember!(Q.cache,key,value,Q.max_cache_entries)
end

# C[w,w,w,0] is known directly from the root data. For d>0 the
# unit equation is x/self_restrictions[w] + remainder = 0.
# Evaluate only the scalar remainder with x=0. Branches unable to reach
# this pivot use the ordinary shared cache instead of a polynomial memo.
function _qs_diagonal(Q, target, d)
  all(iszero,d) && return Q.self_restrictions[target]
  pivot_key = (target,d)
  haskey(Q.diagonal_cache,pivot_key) && return Q.diagonal_cache[pivot_key]
  memo = Dict{Tuple{Int,Int,Int},elem_type(Q.field)}()
  function remainder(u,v,w,dd)
    u == 0 && return zero(Q.field)
    if dd != d || !_qs_below(Q,u,target) || !_qs_below(Q,v,target) || !_qs_below(Q,target,w)
      return _qs_coefficient(Q,u,v,w,dd)
    end
    u == v == w == target && return zero(Q.field)
    u > v && ((u,v) = (v,u))
    key = (u,v,w)
    haskey(memo,key) && return memo[key]
    residual = _qs_recur(Q,u,v,w,dd,remainder)
    return _qs_remember!(memo,key,residual,Q.max_cache_entries)
  end
  value = -_qs_recur(Q,Q.unit,target,target,d,remainder)*Q.self_restrictions[target]
  return _qs_remember!(Q.diagonal_cache,pivot_key,value,Q.max_cache_entries)
end

"""
    quantum_schubert_coefficient(Q, u, v, w, degree=0)

Return the ordinary coefficient of `q^degree sigma_w` in
`sigma_u ⋆ sigma_v` over `Q.field`: an exact rational integer or a residue modulo `Q.prime`.
Indices are context vertex
indices or representative strings. Degrees are nonnegative tuples ordered by
`Q.omitted_roots`; an integer is accepted for Picard rank one (zero for any rank).
"""
function quantum_schubert_coefficient(Q::QuantumSchubertContext,u,v,w,degree=0)
  u,v,w = _qs_vertex(Q,u),_qs_vertex(Q,v),_qs_vertex(Q,w)
  d = _qs_degree(Q,degree)
  Q.codimensions[u]+Q.codimensions[v] == Q.codimensions[w]+_qs_degree_value(Q,d) || return zero(Q.field)
  if Q.codimensions[u] == 1 || Q.codimensions[v] == 1
    a,b = Q.codimensions[u] == 1 ? (u,v) : (v,u)
    i = Int(only(word(Q.representatives[a])))
    return get(quantum_chevalley_product(Q,i,b),(w,d),zero(Q.field))
  end
  c = _qs_coefficient(Q,u,v,w,d)
  isnothing(Q.prime) && !(denominator(c) == 1 && c >= 0) && error("reconstruction produced a nonintegral or negative ordinary coefficient")
  return c
end

"""
    quantum_schubert_product(Q, u, v; max_degree=nothing)

Compute an ordinary quantum basis product as a sparse dictionary
`(output_vertex, degree_tuple) => coefficient` over `Q.field`.
By default all possible degrees are
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
  result = Dict{Tuple{Int,_QSDegree},elem_type(Q.field)}()
  total = Q.codimensions[u]+Q.codimensions[v]
  for d in _qs_degrees(Q,total)
    isnothing(cap) || all(d[i] <= cap[i] for i in eachindex(d)) || continue
    outdegree = total-_qs_degree_value(Q,d)
    0 <= outdegree < length(Q.degree_vertices) || continue
    for w in Q.degree_vertices[outdegree+1]
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

# Classical C[u,w,w,0] = sigma_u|_w. Backward Billey recursion visits
# only right descents of the requested u, rather than all subwords of w
# or all upward Chevalley paths from u. It needs no GKM classes or graph.
function _qs_restriction(Q,u,w)
  u == Q.unit && return one(Q.field)
  u == w && return Q.self_restrictions[w]
  cache_key = (u,w)
  haskey(Q.restriction_cache,cache_key) && return Q.restriction_cache[cache_key]
  target = Q.representatives[u]
  letters = word(Q.representatives[w])
  reflections = reflection.(simple_roots(root_system(parent(target))))
  memo = Dict{Tuple{typeof(target),Int},elem_type(Q.field)}()
  function restriction(x,k)
    len = length(x)
    len == 0 && return one(Q.field)
    len > k && return zero(Q.field)
    key = (x,k)
    haskey(memo,key) && return memo[key]
    residual = restriction(x,k-1)
    xs = x*reflections[Int(letters[k])]
    if length(xs) == len-1
      residual += Q.inversion_values[w][k]*restriction(xs,k-1)
    end
    return _qs_remember!(memo,key,residual,Q.max_cache_entries)
  end
  value = restriction(target,length(letters))
  return _qs_remember!(Q.restriction_cache,cache_key,value,Q.max_cache_entries)
end
