@doc raw"""
    Seidel_element(G::AbstractGKM_graph, w::AbstractAlgebra.Generic.FreeModuleElem{R}, cMax::Int64) where R<:GKM_weight_type

Return the Seidel element on `G` associated to the homomorphism $\iota\colon \mathbb{C}^\times\rightarrow T$ defined by `w`.
This calculates the contribution of curve classes with Chern number at most `cMax`.

!!! warning
    This requires `is_strictly_NEF(G) == true` as otherwise there might be infinitely many contributing curve classes.

# Curve classes
The curve classes of the Seidel space are represented as $(\beta, n)$ where $\beta$ is a curve class of $X$ and the integer $n$
is the degree of the curve projected along $S_X\rightarrow\mathbb{P}^1$.

 - The inclusions of the fibre $X$ over $0$ and $\infty$ simply send $\beta \mapsto (\beta,0)$.
 - Each fixed point $v\in X^T$ gives an invariant curve that is a section of $S_X\rightarrow \mathbb{P}^1$.
    Its curve class is $(\beta_v,1)$, where the $\beta_v$ are only uniquely determined up to an overall constant.
    The Seidel space used in this function is normalized such that $\beta_1=0$, where $1$ is the fixed point represented by the GKM graph vertex with index `1`.

# Example
```jldoctest
julia> P1 = projective_space(GKM_graph, 1);

julia> t1, t2 = gens(P1.M);

julia> cMax = 6;

julia> S1 = Seidel_element(P1, t1, cMax)
(0, -t1 + t2) q^(-1)

julia> S2 = Seidel_element(P1, t2, cMax)
(t1 - t2, 0) q^(0)

julia> S3 = Seidel_element(P1, t1+t2, cMax)
(1, 1) q^(0)

julia> S1 * S2 == S3
true

julia> S4 = Seidel_element(P1, -2*t1, cMax)
(t1^2 - 2*t1*t2 + t2^2, 0) q^(0)
 + (1, 1) q^(1)

julia> S4 * S1 * S1
(1, 1) q^(0)
```
"""
function Seidel_element(G::AbstractGKM_graph, # TODO: calculate maximum chern number of nonzero contribution to get rid of cMax!
    w::AbstractAlgebra.Generic.FreeModuleElem{R},
    cMax::Int64
  ) where R<:GKM_weight_type

  @req is_strictly_nef(G) "G is required to be strictly NEF for the full Seidel element."
  SG = Seidel_space(G, w) # by default, we take the basepoint 1 here.

  nv = n_vertices(G.g)

  # calculate minimum Chern number of v_0 -> v_inf for v a vertex of G:
  cMin = ZZ(0)
  for v in 1:nv
    cMin = min(cMin, chern_number(Edge(v, v+nv), SG))
  end

  #initialize classes to integrate over:
  res = Dict{CurveClass_type, Any}()
  P = [ev(1, point_class(nv+i, SG)) for i in 1:nv]
  coeffRing = G.equivariantCohomology.coeffRingLocalized
  g = gens(coeffRing)
  SGCC = GKM_second_homology(SG)
  GH2proj = SGCC.verticalProjection

  for c in cMin:cMax
    for b in _effectiveSectionClassesWithChernNumber(SG, c)
      # println("Chern number $c, curve class $b -> $(GH2proj(b))")
      GW = gromov_witten(SG, b, 1, P; check_degrees=true, show_bar=false) #TODO: remove this and have global flag for checking
      @req !haskey(res, GH2proj(b)) "Cannot have two distinct section classes restricting to the same vertical curve class!"
      res[GH2proj(b)] = [evaluate(GW[i], vcat(g, [zero(coeffRing)])) for i in 1:nv]
      # println("$(GH2proj(b)) -> $(res[GH2proj(b)])")
    end
  end

  return QH_class(G, res)
end
