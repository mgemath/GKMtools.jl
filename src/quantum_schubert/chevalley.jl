function _qs_minimal(w, levi_reflections)
  while true
    changed = false
    for s in levi_reflections
      ws = w*s
      if length(ws) < length(w)
        w = ws
        changed = true
      end
    end
    changed || return w
  end
end

"""
    quantum_schubert_context(G)
    quantum_schubert_context(R::RootSystem, levi_indices)

Prepare sparse quantum Chevalley data and lazy associativity reconstruction.
The graph overload preserves vertex numbering. The root-system overload
enumerates minimal coset representatives without a GKM graph or the ambient
Weyl group; its ordering is by length and word. `levi_indices` are the simple
roots retained in the parabolic, as in `generalized_gkm_flag`.
"""
function quantum_schubert_context(R::RootSystem, levi_indices)
  levi = sort!(unique(Int.(collect(levi_indices))))
  all(i -> 1 <= i <= rank(R), levi) || throw(ArgumentError("invalid Levi simple root"))
  W = weyl_group(R)
  reflections = reflection.(simple_roots(R))
  lr = reflections[levi]
  reps = [one(W)]
  seen = Set(reps)
  for w in reps
    for s in reflections
      v = _qs_minimal(s*w, lr)
      if v ∉ seen
        push!(seen, v)
        push!(reps, v)
      end
    end
  end
  sort!(reps; by=w -> (length(w), Tuple(Int.(word(w)))))
  return _qs_context(R, levi, reps)
end

function quantum_schubert_context(G::AbstractGKMGraph{R,V,F}) where {R,V<:GeneralizedFlagVertex,F}
  reps = collect(flag.(vertices_structure(G)))
  roots = root_system(parent(first(reps)))
  present = Set(reps)
  levi = [i for i in 1:rank(roots) if reflection(simple_root(roots, i)) ∉ present]
  lr = reflection.(simple_roots(roots))[levi]
  for w in reps, s in reflection.(simple_roots(roots))
    _qs_minimal(s*w, lr) in present || throw(ArgumentError("vertices must form a complete homogeneous G/P"))
  end
  return _qs_context(roots, levi, reps)
end

quantum_schubert_context(::AbstractGKMGraph) = throw(ArgumentError("quantum Schubert calculus requires Weyl-backed homogeneous vertices"))

function _qs_context(R, levi, reps)
  omitted = setdiff(collect(1:rank(R)), levi)
  n = length(reps)
  lengths = length.(reps)
  unit = only(findall(iszero, lengths))
  index = Dict(w => i for (i,w) in enumerate(reps))
  sr = simple_roots(R)
  refl = reflection.(sr)
  lr = refl[levi]
  roots = collect(positive_roots(R))
  root_ids = [k for k in eachindex(roots) if any(i -> !iszero(Oscar.coefficients(roots[k])[i]), omitted)]
  total = [QQ(0) for _ in 1:rank(R)]
  for k in root_ids, i in 1:rank(R)
    total[i] += Oscar.coefficients(roots[k])[i]
  end
  A = cartan_matrix(R)
  nd = [Int(sum(total[j]*A[i,j] for j in 1:rank(R))) for i in omitted]
  all(>(0), nd) || error("nonpositive anticanonical degree")
  outgoing = [_QSEdge[] for _ in reps]
  incoming = [_QSEdge[] for _ in reps]
  zd = Tuple(zeros(Int, length(omitted)))
  for k in root_ids
    alpha = roots[k]
    c = Oscar.coefficients(positive_coroot(R, k))
    d = Tuple(Int(c[i]) for i in omitted)
    coefficients = collect(d)
    height = sum((nd[i]*d[i] for i in eachindex(nd)); init=0)
    s = reflection(alpha)
    for (u,w) in enumerate(reps)
      ws = w*s
      v = get(index, ws, 0)
      if v != 0 && lengths[v] == lengths[u]+1
        push!(outgoing[u], _QSEdge(v, zd, coefficients))
        push!(incoming[v], _QSEdge(u, zd, coefficients))
      end
      v = index[_qs_minimal(ws, lr)]
      if lengths[v] == lengths[u]+1-height
        push!(outgoing[u], _QSEdge(v, d, coefficients))
        push!(incoming[v], _QSEdge(u, d, coefficients))
      end
    end
  end
  # Billey diagonal for D = sum of divisors, in simple-root coordinates.
  # A sufficiently large radix separates all forms and is positive on roots.
  forms = [zeros(BigInt, rank(R)) for _ in reps]
  for (u,w) in enumerate(reps)
    prefix = one(parent(w))
    for letter in word(w)
      i = Int(letter)
      if i in omitted
        beta = Oscar.coefficients(sr[i]*inv(prefix))
        for j in 1:rank(R)
          forms[u][j] += BigInt(ZZ(beta[j]))
        end
      end
      prefix *= refl[i]
    end
  end
  radix = 2*maximum((abs(x) for f in forms for x in f); init=BigInt(0))+2
  diagonal = [QQ(sum((f[i]*radix^(i-1) for i in eachindex(f)); init=BigInt(0))) for f in forms]
  length(Set(diagonal)) == n || error("equivariant divisor specialization did not separate vertices")
  return QuantumSchubertContext(reps, lengths, omitted, nd, unit,
    outgoing, incoming, diagonal, Dict{_QSKey,QQFieldElem}())
end

"""
    quantum_chevalley_product(Q, simple_root_index, u)

Ordinary quantum product of `sigma_u` with the divisor for an omitted simple
root. Returns a sparse dictionary `(vertex, degree_tuple) => QQ`.
"""
function quantum_chevalley_product(Q::QuantumSchubertContext, i::Integer, u)
  u = _qs_vertex(Q, u)
  j = findfirst(==(i), Q.omitted_roots)
  isnothing(j) && throw(ArgumentError("the divisor must correspond to an omitted simple root"))
  result = Dict{Tuple{Int,_QSDegree},QQFieldElem}()
  for e in Q.outgoing[u]
    c = e.coefficients[j]
    iszero(c) && continue
    key = (e.vertex, e.degree)
    result[key] = get(result, key, QQ(0)) + c
  end
  return result
end
