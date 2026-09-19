@doc raw"""
    conjecture_O_eigenvalues(G::AbstractGKMGraph; printData::Bool=true)

Return the eigenvalues of quantum multiplication by ``c_1^T(TX)`` at ``q=1``
and at the non-equivariant limit ``t=0``.

The characteristic polynomial is computed before taking the limit, so this is
well-defined even when the standard-basis matrix contains rational functions
in the equivariant parameters. The graph must be compact and strictly nef.
When `printData=true`, the specialized characteristic polynomial is printed.

# Example
```jldoctest conjecture_O_eigenvalues
julia> P1 = projective_space(GKMGraph, 1);

julia> conjecture_O_eigenvalues(P1)
Characteristic poly of c1(TX)* at q=1, t=0:
x^2 - 4
2-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a1: -2.00000}
```
"""
function conjecture_O_eigenvalues(G::AbstractGKMGraph; printData::Bool=true)
  @req is_compact(G) "the graph must be compact to take the non-equivariant limit"
  chi = characteristic_polynomial(c1_at_q1(G))
  chi0 = polynomial(QQ, [0])
  z = zeros(Int, rank_torus(G))
  for i in 0:(length(chi) - 1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  printData && println("Characteristic poly of c1(TX)* at q=1, t=0:\n$chi0")
  roots(QQBar, chi0)
end
