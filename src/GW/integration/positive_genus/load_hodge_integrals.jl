export load_H

function load_H()
  return _load_Hodge_integrals(joinpath(@__DIR__, "../../../../M2"))
end

function load_H(gMax::Int64, nMax::Int64)
  return _load_Hodge_integrals(joinpath(@__DIR__, "../../../../M2"); gMax=gMax, nMax=nMax)
end

# What was the envisaged purpose of the function below?
# function load_H(::Int64)
#   return Dict{HodgeKey, QQFieldElem}()
# end
const _HODGE_INTEGRAL_CACHE = Dict{Tuple{Int,Int},Dict{HodgeKey,QQFieldElem}}()

function _cached_hodge_integrals(max_genus::Int, max_marks::Int)
  return get!(_HODGE_INTEGRAL_CACHE, (max_genus, max_marks)) do
    load_H(max_genus, max_marks)
  end
end
