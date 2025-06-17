struct TreeIt
  n_vert::Int64
end

Base.eltype(::Type{TreeIt}) = Vector{Int64}

Base.@propagate_inbounds function Base.length(TI::TreeIt)::Int64
  n = TI.n_vert
  n < 4 && return 1
  a = A000081(n)
  len = a[n] - sum(i -> a[i] * a[n - i], 1:(n >>> 1))
  if iseven(n)
    a_half = a[n >>> 1]
    len += ((a_half + 1) * a_half) >>> 1
  end
  return len
end

function Base.iterate(TI::TreeIt, state::Tuple=())::Union{Nothing,Tuple{Vector{Int64},Tuple}}
  if isempty(state)
    n = TI.n_vert
    if n < 4
      return map(i -> min(2, i), 1:n), (0,)
    end
    r = n >>> 1 + 1
    l = [1:r; 2:n - r + 1]
    w = [0:r - 1; 1:n - 1]
    return l, (n - 1, l, w, n, n, r, n, iseven(n) ? Float64(n + 1) : Inf,
               r, false, false, false, false)
  end

  q, l, w, n, p, h1, h2, c, r, fixit, needr, needc, needh2 = state
  q == 0 && return nothing

  if c == n + 1 || (p == h2 && ((l[h1] == l[h2] + 1 && n - h2 > r - h1) ||
                                (l[h1] == l[h2] && n - h2 + 1 < r - h1)))
    if l[r] > 3
      p, q = r, w[r]
      h1 == r && (h1 -= 1)
      fixit = true
    else
      p, r, q = r, r - 1, 2
    end
  end

  p <= h1 && (h1 = p - 1)
  if p <= r
    needr = true
  elseif p <= h2
    needh2 = true
  elseif l[h2] == l[h1] - 1 && n - h2 == r - h1
    p <= c && (needc = true)
  else
    c = Inf
  end

  oldp, δ, oldlq, oldwq = p, q - p, l[q], w[q]
  p = Inf

  for i in oldp:n
    l[i] = l[i + δ]
    if l[i] == 2
      w[i] = 1
    else
      p = i
      q = l[i] == oldlq ? oldwq : w[i + δ] - δ
      w[i] = q
    end

    if needr && l[i] == 2
      needr, needh2, r = false, true, i - 1
    end

    if needh2 && l[i] <= l[i - 1] && i > r + 1
      needh2, h2 = false, i - 1
      if l[h2] == l[h1] - 1 && n - h2 == r - h1
        needc = true
      else
        c = Inf
      end
    end

    if needc
      c = l[i] != l[h1 - h2 + i] - 1 ? Inf : Float64(i + 1)
      needc = false
    end
  end

  if fixit
    r = n - h1 + 1
    for i in r + 1:n
      l[i] = i - r + 1
      w[i] = i - 1
    end
    w[r + 1] = 1
    h2, p, q, c = n, n, n - 1, Inf
  elseif p == Inf
    p = l[oldp - 1] != 2 ? oldp - 1 : oldp - 2
    q = w[p]
    if needh2
      h2 = n
      c = (l[h2] == l[h1] - 1 && h1 == r) ? Float64(n + 1) : Inf
    end
  end

  return l, (q, l, w, n, p, h1, h2, c, r, false, false, false, false)
end

function LStoGraph(ls::Vector{Int64})::Graph{Undirected}
  n = length(ls)
  g = Graph{Undirected}(n)
  for v in 2:n
    p = findlast(i -> i < v && ls[i] == ls[v] - 1, eachindex(ls))
    add_edge!(g, p, v)
  end
  return g
end

function A000081(n::Int64)::Vector{Int64}
  n < 3 && return ones(Int64, n)
  a = Vector{Int64}(undef, n)
  a[1:4] = (1, 1, 2, 4)
  for i in 5:n
    acc = a[i - 1]
    for k in 2:i - 1
      sumdiv = k * a[k] + 1
      for d in 2:k - 1
        k % d == 0 && (sumdiv += d * a[d])
      end
      acc += sumdiv * a[i - k]
    end
    a[i] = acc ÷ (i - 1)
  end
  return a
end

function A000055(n::Int64)::Vector{Int64}
  n < 4 && return ones(Int64, n)
  a = A000081(n)
  len = Vector{Int64}(undef, n)
  len[1] = 1
  for j in 2:n
    len[j] = a[j] - sum(i -> a[i] * a[j - i], 1:(j >>> 1))
  end
  for j in 2:2:n
    a_half = a[j >>> 1]
    len[j] += ((a_half + 1) * a_half) >>> 1
  end
  return len
end
