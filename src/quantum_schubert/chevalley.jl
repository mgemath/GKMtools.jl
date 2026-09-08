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
    quantum_schubert_context(G; p=nothing, max_cache_entries=200_000)
    quantum_schubert_context(R::RootSystem, levi_indices; p=nothing, max_cache_entries=200_000)

Prepare sparse quantum Chevalley data and lazy associativity reconstruction.
The graph overload preserves vertex numbering. The root-system overload
enumerates minimal coset representatives without a GKM graph or the ambient
Weyl group; its ordering is by length and word. `levi_indices` are the simple
roots retained in the parabolic, as in `generalized_gkm_flag`.

Omit `p` for exact rational computation, or supply a prime to reconstruct
directly over GF(p). The prime must exceed the basis size and admit a
separating equivariant specialization; invalid choices throw ArgumentError.
Each memo is bounded by `max_cache_entries` (zero disables memoization).
This limits retained cache entries, not total RAM.
"""
function quantum_schubert_context(R::RootSystem, levi_indices; p=nothing, max_cache_entries::Int=200_000)
  _qs_field(p) # validate before enumerating representatives
  max_cache_entries >= 0 || throw(ArgumentError("max_cache_entries must be nonnegative"))
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
  return _qs_context(R, levi, reps; p, max_cache_entries)
end

function quantum_schubert_context(G::AbstractGKMGraph{R,V,F}; p=nothing, max_cache_entries::Int=200_000) where {R,V<:GeneralizedFlagVertex,F}
  _qs_field(p)
  max_cache_entries >= 0 || throw(ArgumentError("max_cache_entries must be nonnegative"))
  reps = collect(flag.(vertices_structure(G)))
  roots = root_system(parent(first(reps)))
  present = Set(reps)
  levi = [i for i in 1:rank(roots) if reflection(simple_root(roots, i)) ∉ present]
  lr = reflection.(simple_roots(roots))[levi]
  for w in reps, s in reflection.(simple_roots(roots))
    _qs_minimal(s*w, lr) in present || throw(ArgumentError("vertices must form a complete homogeneous G/P"))
  end
  return _qs_context(roots, levi, reps; p, max_cache_entries)
end

quantum_schubert_context(::AbstractGKMGraph; kwargs...) = throw(ArgumentError("quantum Schubert calculus requires Weyl-backed homogeneous vertices"))

function _qs_context(R, levi, reps; p=nothing, max_cache_entries::Int=200_000)
  omitted = setdiff(collect(1:rank(R)), levi)
  n = length(reps)
  isnothing(p) || p > n || throw(ArgumentError("p must exceed the number of Schubert classes for a separating specialization"))
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
  K = _qs_field(p)
  values, diagonal = _qs_specialization(K, p, forms, roots)
  # Classical C[w,w,w,0] is Billey's diagonal restriction: the product
  # of the inversion roots in a reduced word. This also supplies the
  # inverse path-sum coefficient in positive-degree diagonal reconstruction.
  restrictions = [one(K) for _ in reps]
  inversions = [elem_type(K)[] for _ in reps]
  for (u,w) in enumerate(reps)
    prefix = one(parent(w))
    for letter in word(w)
      i = Int(letter)
      beta = Oscar.coefficients(sr[i]*inv(prefix))
      value = sum((K(ZZ(beta[j]))*values[j] for j in eachindex(values)); init=zero(K))
      push!(inversions[u],value)
      restrictions[u] *= value
      prefix *= refl[i]
    end
  end
  blocks = [findall(==(d),lengths) for d in 0:maximum(lengths)]
  # Classical Chevalley covers generate Bruhat order. Bitsets avoid
  # repeatedly comparing reduced words in the innermost recurrence.
  # Limit this optional index to 64 MiB for much larger spaces.
  below = if BigInt(n)*n <= 8*64*1024^2
    intervals = [falses(n) for _ in reps]
    for u in sortperm(lengths)
      intervals[u][u] = true
      for e in outgoing[u]
        all(iszero,e.degree) || continue
        intervals[e.vertex] .|= intervals[u]
      end
    end
    intervals
  else
    nothing
  end
  T = elem_type(K)
  return QuantumSchubertContext(reps, lengths, omitted, nd, unit,
    outgoing, incoming, K, isnothing(p) ? nothing : BigInt(p), diagonal,
    restrictions, inversions, Dict{Tuple{Int,Int},T}(), Dict{_QSKey,T}(), Dict{Tuple{Int,_QSDegree},T}(),
    max_cache_entries, below, blocks)
end

"""
    quantum_chevalley_product(Q, simple_root_index, u)

Ordinary quantum product of `sigma_u` with the divisor for an omitted simple
root. Returns a sparse dictionary `(vertex, degree_tuple) => coefficient`
over `Q.field` (QQ or GF(p)).
"""
function quantum_chevalley_product(Q::QuantumSchubertContext, i::Integer, u)
  u = _qs_vertex(Q, u)
  j = findfirst(==(i), Q.omitted_roots)
  isnothing(j) && throw(ArgumentError("the divisor must correspond to an omitted simple root"))
  result = Dict{Tuple{Int,_QSDegree},elem_type(Q.field)}()
  for e in Q.outgoing[u]
    c = e.coefficients[j]
    iszero(c) && continue
    key = (e.vertex, e.degree)
    result[key] = get(result, key, zero(Q.field)) + Q.field(c)
  end
  filter!(pair -> !iszero(last(pair)), result)
  return result
end

# Exact computation uses an injective radix. Modulo p, try deterministic
# specializations and certify separation and nonvanishing of every root.
# Failure is explicit: division by zero must never silently corrupt output.
function _qs_specialization(K, p, forms, roots)
  radix = 2*maximum((abs(x) for f in forms for x in f); init=BigInt(0))+2
  r = length(first(forms))
  seed = BigInt(1729)
  for attempt in 1:(isnothing(p) ? 1 : 128)
    values = if attempt == 1
      [K(radix^(i-1)) for i in 1:r]
    else
      [begin
        seed = mod(6364136223846793005*seed+1442695040888963407,BigInt(2)^64)
        K(mod(seed,p))
      end for _ in 1:r]
    end
    diagonal = [sum((K(f[i])*values[i] for i in 1:r);init=zero(K)) for f in forms]
    length(Set(diagonal)) == length(forms) || continue
    all(alpha -> !iszero(sum((K(ZZ(Oscar.coefficients(alpha)[i]))*values[i] for i in 1:r);init=zero(K))),roots) || continue
    return values, diagonal
  end
  throw(ArgumentError("no nonsingular equivariant specialization found modulo p; choose a larger prime"))
end
