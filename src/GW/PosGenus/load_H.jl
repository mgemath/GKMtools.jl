export load_H

function load_H()
  return _load_Hodge_integrals(joinpath(@__DIR__ , "../../../M2/"))
end

function load_H(::Int64)
  return Dict{HodgeKey, QQFieldElem}()
end