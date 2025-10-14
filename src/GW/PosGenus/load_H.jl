export load_H

function load_H()
  return _load_Hodge_integrals(joinpath(@__DIR__ , "../../../M2/"))
end

function load_H(gMax::Int64, nMax::Int64)
  return _load_Hodge_integrals(joinpath(@__DIR__ , "../../../M2/"); gMax=gMax, nMax=nMax)
end

# What was the envisaged purpose of the function below?
# function load_H(::Int64)
#   return Dict{HodgeKey, QQFieldElem}()
# end