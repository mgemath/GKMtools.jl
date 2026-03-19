@doc raw"""
    Psi(a) -> EquivariantClass

For each index $i$ such that $1\le i \le n$, there is a line bundle on $\overline{\mathcal{M}}_{g,n}(X,\beta)$ such that the fiber at a moduli point is the cotangent bundle of the curve at the $i^\mathrm{th}$ marked point. 
We denote by $\psi_i$ the first Chern class of this line bundle. In order to compute invariants involving ${\psi_1}^{a_1}\cdots {\psi_n}^{a_n}$, for some nonnegative integers  $a_1,\ldots, a_n$, we write `Psi(a_1,...,a_n)`.


# Examples
Let $G$ be the GKM graph of the Hirzebruch surface $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(1))$, let $\beta$ the class of the fiber. The invariant

```math
\int_{\overline{\mathcal{M}}_{0,2}(G, \beta)}\mathrm{ev}_{1}^{*}([\mathrm{pt}])\cdot\psi_{1}^{0}\psi_{2} = -1,
```

can be computed as following.
```jldoctest
julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 1));

julia> P = ev(1, point_class(G, 1)) * Psi(0,1);

julia> beta = curve_class(G, "1", "4"); # beta is a fiber of the map G -> P^1

julia> gromov_witten(G, beta, 2, P; show_bar=false)
-1
```

Let us give an example in positive genus. Let `P1` be the GKM graph of $\mathbb{P}^1$, and let $\beta$ be the class of the line. Choose $n$ and 
$g$ such that $d=n+1-g>0$. The invariant 

```math
\int_{\overline{\mathcal{M}}_{g,n}(\mathbb{P}^1, d\beta)}\prod_{i=1}^n \left(\mathrm{ev}_{i}^{*}([\mathrm{pt}]) \psi_{i}^2 \right),
```
can be computed as following:
```jldoctest Psi 2
julia> P1 = projective_space(GKM_graph, 1); 

julia> beta = curve_class(P1, "1", "2"); # beta is the class of the line

julia> w = point_class(1, P1);

julia> n = 4; g = 2; # so that d > 0

julia> d = n + 1 - g;

julia> P = prod(i -> ev(i, w), 1:n); # = ev(1, w) * ev(2, w) * ev(3, w) * ev(4, w)
```
Now we multiply `P` by the psi-classes:
```jldoctest Psi 2
julia> P *= Psi([2 for _ in 1:n]); # = Psi(2,2) * ev(1, w)*ev(2, w)*ev(3, w)*ev(4, w)

julia> gromov_witten(P1, d*beta, n, P; g = g, show_bar = false)
263//96
```
Results match those in [MR4028099](@cite).

To see a non-compact example, let $X:=\text{Tot}(\mathcal{O}_{\mathbb{P}^2}(-1)\oplus\mathcal{O}_{\mathbb{P}^2}(-2))$ and let us compute
```math
\int_{\overline{\mathcal{M}}_{0,2}(X, d\beta)}\psi_1\psi_2^2
```
for $1\le d\le 4$, where $\beta$ is the hyperplane class.

```jldoctest
julia> G = total_space(vector_bundle_O(2, [-1, -2]));

julia> beta = curve_class(G, Edge(1, 2));

julia> P = Psi(1,2);

julia> for d in 1:4
           println(gromov_witten(G, d*beta, 2, P; show_bar=false))
       end
0
1
-3//2
401//144
```

Finally, let $X:=\text{Tot}(\mathcal{O}_{\mathbb{P}^1}(-1)\oplus\mathcal{O}_{\mathbb{P}^1}(-1))$ and let us compute
```math
\int_{\overline{\mathcal{M}}_{g,2}(X, d\beta)}\psi_1\psi_2
```
for $1\le d\le 3$ and $0\le g\le 3$, where $\beta$ is the hyperplane class.

```jldoctest
julia> G = total_space(vector_bundle_O(1, [-1, -1]));

julia> beta = curve_class(G, Edge(1, 2));

julia> P = Psi(1, 1);

julia> for g in 0:3
           for d in 1:3
               print("g=$g, d=$d: ")
               println(gromov_witten(G, d*beta, 2, P; g=g, show_bar=false))
           end
       end
g=0, d=1: 2
g=0, d=2: 1//4
g=0, d=3: 2//27
g=1, d=1: 0
g=1, d=2: 0
g=1, d=3: 0
g=2, d=1: 1//40
g=2, d=2: 1//20
g=2, d=3: 3//40
g=3, d=1: 5//1512
g=3, d=2: 5//189
g=3, d=3: 5//56
```
"""
function Psi(a::Int64)::EquivariantClass
  # The following threw a type error.
  # rule = :(_Psi(dt, $a))
  # return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, [Int64(a)], Int64[])
  return Psi([a])
end

function Psi(a::Vector{Int64})::EquivariantClass
  rule = :(_Psi(dt, $a))
  return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, a, Int64[])
end

function Psi(a::Int...)::EquivariantClass
  # rule = :(_Psi(dt, $a))
  # return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, collect(a), Int64[])
  return Psi(collect(a))
end

function _Psi(dt::GW_decorated_tree, a::Int64...)
  return _Psi(dt, [a])
end

function _Psi(dt::GW_decorated_tree, a::Tuple{Vararg{Int64}})
  return _Psi(dt, [a...])
end

function _Psi(dt::GW_decorated_tree, a::Vector{Int64})
  findfirst(x -> x > 0, a) === nothing && return 1 # F(1) #if all of them are zero or a is empty
  g = dt.tree
  marks = dt.marks

  ans = 1 #F(1)

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

    n > 2 && Sum_ai > n - 3 && return 0# F(0)

    #If no previous condition holds, then n>1
    if n == 2 #necessary |S_v| == 1
      M = (-1)^a_v[1]
    else # n>2 and Sum_ai <= n - 3
      M = multinomial(n - 3 - Sum_ai, a_v...)
    end

    local s1 = 0# F(0)

    for w in all_neighbors(g, v)
      e = Edge(v, w)
      # wev = weight_class(imageOf(e, dt), R) // edgeMult(e, dt)
      s1 += edgeMult(e, dt)//weight_class(imageOf(e, dt), dt.gkm) #  1 // wev # flag-compatible since only concerns edges of the decorated tree.
    end
    ans *= M * (s1^(-Sum_ai))
  end

  return ans
end

function _Psi(dg::GW_decorated_graph, a)
  return 1
end