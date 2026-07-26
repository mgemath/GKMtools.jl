function Psi(a::Int64)::EquivariantClass
  # The following threw a type error.
  # rule = :(_Psi(dt, $a))
  # return EquivariantClass(rule, (dt -> _Psi(dt, a)), false, false, [Int64(a)], Int64[])
  return Psi([a])
end

function Psi(a::Vector{Int64})::EquivariantClass
  rule = :(_Psi(dt, $a))
  return EquivariantClass(rule, dt -> _Psi(dt, a), false, false, a, Int64[])
end

function Psi(a::Int...)::EquivariantClass
  # rule = :(_Psi(dt, $a))
  # return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, collect(a), Int64[])
  return Psi(collect(a))
end

function _Psi(dt, a::Int64...)
  return _Psi(dt, [a])
end

function _Psi(dt, a::Tuple{Vararg{Int64}})
  return _Psi(dt, [a...])
end

function _Psi(dt, a::Vector{Int64})
  t = dt.context.t
  scalar_one = one(t[1]) // one(t[1])
  hasproperty(dt, :genus) && return scalar_one
  all(iszero, a) && return scalar_one
  g = dt.tree
  marks = dt.marks

  ans = scalar_one

  # local q1::fmpq = fmpq(1)
  # local temp1::fmpq = fmpq(1)
  local Sum_ai::Int64
  local n::Int64
  local M::Int64
  # local d = Dict(edges(g) .=> weights) #assign weights to edges
  local inv_marks::Dict{Int64,Vector{Int64}} = invert_marks(marks, nv(g))

  for v in 1:nv(g)
    a_v = Int64[]
    for i in inv_marks[v]
      (i > length(a) || a[i] == 0) && continue
      push!(a_v, a[i])
    end

    Sum_ai = sum(a_v)
    Sum_ai == 0 && continue #if S contains only zeros, or it is empty, continue

    n = length(all_neighbors(g, v)) + length(inv_marks[v])

    n > 2 && Sum_ai > n - 3 && return zero(scalar_one)

    #If no previous condition holds, then n>1
    if n == 2 #necessary |S_v| == 1
      M = (-1)^a_v[1]
    else # n>2 and Sum_ai <= n - 3
      M = multinomial(n - 3 - Sum_ai, a_v...)
    end

    s1 = zero(scalar_one)

    for w in all_neighbors(g, v)
      e = Edge(v, w)
      # wev = weight_class(imageOf(e, dt), R) // edgeMult(e, dt)
      s1 += edgeMult(e, dt) // _weight_class(dt.gkm, imageOf(e, dt), t) #  1 // wev # flag-compatible since only concerns edges of the decorated tree.
    end
    ans *= M * (s1^(-Sum_ai))
  end

  return ans
end

