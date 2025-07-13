@doc raw"""
  Psi(a) -> EquivariantClass

For each index $i$ such that $0\le i \le n$, there is a line bundle on $\overline{\mathcal{M}_{0,n}}(X,\beta)$ such that the fiber at a moduli point is the cotangent bundle of the curve at the $i^\mathrm{th}$ marked point. 
We denote by $\psi_i$ the first Chern class of this line bundle. In order to compute invariants involving ${\psi_1}^{a_1}\cdots {\psi_n}^{a_n}$, for some nonnegative integers  $a_1,\ldots, a_n$, we write `Psi(a_1,...,a_n)`

# Example
Let $G$ be the GKM graph of the Hirzebruch surface $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(1))$, let $\beta$ the class of the fiber. The invariant

```math
\int_{\overline{M}_{0,2}(G, \beta)}\mathrm{ev}_{1}^{*}([\mathrm{pt}])\cdot\psi_{1}^{0}\psi_{2} = -1,
```

can be computed as following."""
# ```jldoctest ev
# julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 1));

# julia> P = ev(1, point_class(G, 1)) * Psi(0,1);

# julia> beta = curve_class(G, "1", "4"); # beta is a fiber of the map G -> P^1

# julia> gromov_witten(G, beta, 2, P; show_bar=false)
# -1
# ```
# """
function Psi(a)::EquivariantClass

  rule = :(_Psi(dt, $a))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function Psi(a::Int...)::EquivariantClass

  rule = :(_Psi(dt, $a))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _Psi(dt::GW_decorated_tree, a::Int64...)
  return _Psi(dt, [a])
end

function _Psi(dt::GW_decorated_tree, a::Tuple{Vararg{Int64}})
  return _Psi(dt, [a...])
end

function _Psi(dt::GW_decorated_tree, a::Vector{Int64})

  ONE = one(dt.gkm.equivariantCohomology.coeffRing) # F(1)
  ZERO = zero(dt.gkm.equivariantCohomology.coeffRing) # F(0)

  findfirst(x -> x > 0, a) === nothing && return ONE # F(1) #if all of them are zero or a is empty
  g = dt.tree
  marks = dt.marks


  ans::QQMPolyRingElem = ONE #F(1)

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

    n > 2 && Sum_ai > n - 3 && return ZERO # F(0)

    #If no previous condition holds, then n>1
    if n == 2 #necessary |S_v| == 1
      M = (-1)^a_v[1]
    else # n>2 and Sum_ai <= n - 3
      M = multinomial((n - 3 - Sum_ai, a_v...,))
    end


    local s1 = ZERO# F(0)

    for w in all_neighbors(g, v)
      e = Edge(v,w)
      # wev = weight_class(imageOf(e, dt), R) // edgeMult(e, dt)
      s1 += edgeMult(e, dt) // weight_class(imageOf(e, dt), dt.gkm) #  1 // wev 
    end
    ans *= M * (s1^(-Sum_ai))
  end

  return ans
end
