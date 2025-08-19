import Base: *, //, +, -, one, zero, ==
import Oscar: is_zero

_CohomType = Union{AbstractAlgebra.Generic.FreeModuleElem, FreeModElem, Array}

@doc raw"""
    QH_class(G::AbstractGKM_graph, class; beta::Union{Nothing, CurveClass_type} = nothing)

Turn the given equivariant cohomology class `class` into an equivariant quantum cohomology class on `G`.
The optional argument `beta` can be used to multiply the result by the coefficient $q^\beta$ for a curve class $\beta$.

# Example
```jldoctest QH_class
julia> P2 = projective_space(GKM_graph, 2);

julia> QH_class(P2, point_class(P2, 1))
(t1^2 - t1*t2 - t1*t3 + t2*t3, 0, 0) q^(0)

julia> QH_class(P2, point_class(P2, 1); beta = curve_class(P2, Edge(1, 2)))
(t1^2 - t1*t2 - t1*t3 + t2*t3, 0, 0) q^(1)

julia> (t1, t2, t3) = gens(P2.equivariantCohomology.coeffRing); # hyperplane class

julia> QH_class(P2, [t1, t2, t3])
(t1, t2, t3) q^(0)
```
"""
function QH_class(G::AbstractGKM_graph, class; beta::Union{Nothing, CurveClass_type} = nothing)
  H2 = GKM_second_homology(G)
  if isnothing(beta)
    beta = zero(H2.H2)
  end
  return QH_class(G, Dict{CurveClass_type, Any}([beta => class]))
end

function QH_class(G::AbstractGKM_graph, classes::Dict{CurveClass_type, Any})
  H2 = GKM_second_homology(G)
  nv = n_vertices(G.g)
  coeffs = Dict{CurveClass_type, QH_coeff_type}()
  Rloc = G.equivariantCohomology.coeffRingLocalized
  g = gens(G.equivariantCohomology.cohomRingLocalized)
  for beta in keys(classes)
    @req parent(beta) == H2.H2 "Beta belongs to the wrong second homology group."
    c0 = zero(G.equivariantCohomology.cohomRingLocalized)
    for i in 1:nv
      c0 += Rloc(classes[beta][i]) * g[i]
    end
    coeffs[beta] = c0
  end
  return _QH_remove_zero_coeffs!(QHRingElem(G, coeffs))
end

function _QH_remove_zero_coeffs!(c::QHRingElem)::QHRingElem
  for b in keys(c.coeffs)
    if is_zero(c.coeffs[b])
      delete!(c.coeffs, b)
    end
  end
  return c
end

function +(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c1.coeffs)
    res.coeffs[b] = c1.coeffs[b]
  end
  for b in keys(c2.coeffs)
    if haskey(res.coeffs, b)
      res.coeffs[b] += c2.coeffs[b]
    else
      res.coeffs[b] = c2.coeffs[b]
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function +(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH + c2
end

function +(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 + c2QH
end

function +(c1::QHRingElem, c2)::QHRingElem
  return c1 + one(c2) * c2
end

function +(c1, c2::QHRingElem)::QHRingElem
  return one(c1)*c1 + c2
end

# Here, s should be a scalar, i.e. anything that embeds into the localized coefficient ring H_T(*)_loc.
function *(s, c::QHRingElem)::QHRingElem
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c.coeffs)
    res.coeffs[b] = (s//1) * c.coeffs[b]
  end
  return _QH_remove_zero_coeffs!(res)
end

function *(c::QHRingElem, s)::QHRingElem
  return *(s, c)
end

@doc raw"""
    *(c1::QHRingElem, c2::QHRingElem) -> QHRingElem

Multiply the classes `c1` and `c2` using the equivariant quantum product in $QH_T^*(X)$.

!!! warning
    This requires `is_strictly_nef(G)==true` for the underlying GKM graph `G`.
    If this does not hold, there could potentially be infinitely many $\beta$ contributing a non-zero $q^\beta$-term to the quantum product.
    In this case, use `quantum_product` to calculate the coefficient of $q^\beta$ in the quantum product for a specified choice of $\beta$.

# Example
```jldoctest QH_product_*
julia> P2 = projective_space(GKM_graph, 2);

julia> (t1, t2, t3) = gens(P2.equivariantCohomology.coeffRing);

julia> H = QH_class(P2, [t1, t2, t3]) # The equivariant hyperplane class as element of QH_T(X)
(t1, t2, t3) q^(0)

julia> (H - t1) * (H - t2) * (H - t3)
(1, 1, 1) q^(1)

julia> p = QH_class(P2, point_class(P2, 1))
(t1^2 - t1*t2 - t1*t3 + t2*t3, 0, 0) q^(0)

julia> p * H
(t1^3 - t1^2*t2 - t1^2*t3 + t1*t2*t3, 0, 0) q^(0)
 + (1, 1, 1) q^(1)
```
"""
function *(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  #calculate QH structure constants in all beta. This throws an error if G is not strictly NEF.
  SC = QH_structure_constants(c1.gkm; show_progress=false)
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b1 in keys(c1.coeffs)
    for b2 in keys(c2.coeffs)
      for b in keys(SC)
        tmp = quantum_product(c1.gkm, b, c1.coeffs[b1], c2.coeffs[b2])
        k = b1+b2+b
        if haskey(res.coeffs, k)
          res.coeffs[k] += tmp
        else
          res.coeffs[k] = tmp
        end
      end
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function *(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH * c2
end

function *(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 * c2QH
end

function -(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c1.coeffs)
    res.coeffs[b] = c1.coeffs[b]
  end
  for b in keys(c2.coeffs)
    if haskey(res.coeffs, b)
      res.coeffs[b] -= c2.coeffs[b]
    else
      res.coeffs[b] = -c2.coeffs[b]
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function -(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH - c2
end

function -(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 - c2QH
end

function -(c1::QHRingElem, c2)::QHRingElem
  return c1 - one(c1) * c2
end

function -(c1, c2::QHRingElem)::QHRingElem
  return one(c2)*c1 - c2
end

function //(c::QHRingElem, s)::QHRingElem
  @req s != 0 "Cannot divide QHRingElem by zero!"
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c.coeffs)
    res.coeffs[b] = (1//s) * c.coeffs[b]
  end
  return res
end

function one(c::QHRingElem)::QHRingElem
  res = zero(c)
  for g in gens(c.gkm.equivariantCohomology.cohomRingLocalized)
    res += g
  end
  return res
end

function zero(c::QHRingElem)::QHRingElem
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  res.coeffs[zero(GKM_second_homology(c.gkm).H2)] = zero(c.gkm.equivariantCohomology.cohomRingLocalized)
  return res
end

function is_zero(c::QHRingElem)::Bool
  for b in keys(c.coeffs)
    if !is_zero(c.coeffs[b])
      return false
    end
  end
  return true
end

function ==(c1::QHRingElem, c2::QHRingElem)::Bool
  c1.gkm != c2.gkm && return false
  for b in keys(c1.coeffs)
    if haskey(c2.coeffs, b)
      c1.coeffs[b] != c2.coeffs[b] && return false
    else
      !is_zero(c1.coeffs[b]) && return false
    end
  end
  for b in keys(c2.coeffs)
    if haskey(c1.coeffs, b)
      c2.coeffs[b] != c1.coeffs[b] && return false
    else
      !is_zero(c2.coeffs[b]) && return false
    end
  end
  return true
end

function ==(c1::QHRingElem, c2::_CohomType)::Bool
  c2QH = QH_class(c1.gkm, c2)
  return c1 == c2QH
end

function ==(c1::_CohomType, c2::QHRingElem)::Bool
  return c2 == c1
end

function ==(c1::QHRingElem, c2)::Bool
  return c1 == one(c1)*c2
end

function ==(c1, c2::QHRingElem)::Bool
  return c2 == c1
end

function Base.show(io::IO, c::QHRingElem)

  if isempty(c.coeffs)
    print(io, "0")
    return
  end

  l = length(keys(c.coeffs))
  counter = 0
  # no nested printing
  for b in keys(c.coeffs)
    counter += 1
    if isnothing(c.gkm.QH_preferred_basis)
      print(io, c.coeffs[b])
    else
      p = c.gkm.QH_preferred_basis * [c.coeffs[b][i] for i in 1:n_vertices(c.gkm.g)]
      print(io, "(")
      for i in 1:length(p)
        (i > 1) && print(io, ", ")
        if p[i] == 0
          print(io, "0")
        elseif denominator(p[i]) == 1
          print(io, factor(numerator(p[i])))
        else
          print(io, "(")
          print(io, factor(numerator(p[i])))
          print(io, ")//(")
          print(io, factor(denominator(p[i])))
          print(io, ")")
        end
      end
      print(io, ")")
    end
    print(io, " q^")
    print(io, b)
    if counter != l
      if Oscar.is_terse(io)
        print(io, " + ")
      else
        print(io, "\n + ")
      end
    end
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", c::QHRingElem)
  show(io, c)
end