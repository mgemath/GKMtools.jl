# Warning: This is not actually the inverse of the Euler class, as the h classes will be multiplied later.
# Note: this throws an error if the decorated tree is a single vertex with less than 3 vertices.
# But GW_decorated_tree() checks this during construction, so no check is done here.
function Euler_inv(dt::GW_decorated_tree, t::Vector{T}, edge_weight_dict::Dict{Edge, T}, point_weight_dict::Vector{Union{Nothing, T}}; check_degree::Bool=false) where T<:RingElem

  res = one(t[1])
  oldDeg = 0

  for v in 1:n_vertices(dt.tree)

    valv = degree(dt.tree, v)
    e = euler_class(imageOf(v, dt), dt.gkm.equivariantCohomology, t, edge_weight_dict, point_weight_dict)
    #println("e = $e, val = $valv")
    if valv >= 1
      res = res * e^(valv - 1)
    else
      res = res // e
    end

    tmpSum = zero(t[1])

    for v2 in all_neighbors(dt.tree, v)
      e = Edge(v,v2)
      wev = weight_class(imageOf(e, dt), dt.gkm, t, edge_weight_dict) // edgeMult(e, dt)
      res = res // wev
      tmpSum = tmpSum + 1//wev
    end

    e =  valv - 3 +  count(i -> i==v, dt.marks)
    if e >= 0
      res = res * tmpSum^e
    else 
      res = res // (tmpSum^(-e))
    end

    if check_degree
      r = valency(dt.gkm)
      Sv = count(i -> i==v, dt.marks)
      corDeg = r*(valv - 1) + 3 - Sv - 2*valv
      if _get_degree(res) - oldDeg != corDeg
        println("Wrong Euler class degree:")
        println("dt: $dt")
        println("r: $r, v=val(v) = $valv")
        @req false "Euler class has wrong degree!"
      end
      oldDeg += corDeg
    end
  end

  return res
end

# This returns the extra factor for Euler_inv in the fiber direction.
function _Euler_inv_VB(dt::GW_decorated_tree, V::GKM_vector_bundle)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  C = dt.gkm.equivariantCohomology.coeffRing
  res = C(1)//C(1)

  for v in 1:n_vertices(dt.tree)

    valv = degree(dt.tree, v)
    e = _fiber_normal_weight(imageOf(v, dt), V)
    #println("e = $e, val = $valv")
    if valv >= 1
      res = res * e^(valv - 1)
    else
      res = res // e
    end
  end

  return res
end

# Calculate h(epsilon, d) as in [Liu--Sheshmani, Lemma 4.5, p. 16].
function _h(e::Edge, d::Int, con::GKM_connection, R::GKM_cohomology_ring, t::Vector{T}, edge_weight_dict::Dict{Edge, T}; check::Bool=true, check_degrees::Bool=false) where T<:RingElem

  gkm = con.gkm

  if check
    @req con.gkm == R.gkm "GKM connection and cohomology ring don't belong to the same GKM graph"
    @req has_edge(gkm.g, e) "edge not found in GKM graph"
    @req d>0 "d is non-positive"
  end

  we = weight_class(e, R.gkm, t, edge_weight_dict) # weight of the edge e

  res = ( one(t[1]) * ((-1)^d) * ZZ(d)^(2d) ) // ( factorial(ZZ(d))^2 )
  res = res // ( (we)^(2d) )

  for v in all_neighbors(gkm.g, src(e))

    v == dst(e) && continue

    ei = Edge(src(e), v)
    wei = weight_class(ei, R.gkm, t, edge_weight_dict)
    ai = con.a[(e, ei)]
    bFactor = _b(1//d * we, wei, d*ai)
    if check_degrees
      corDeg = -d*ai - 1
      actDeg = _get_degree(bFactor)
      if actDeg != corDeg
        println("b factor has wrong deg for d=$d, ai=$ai, we=$we, wei=$wei.")
        println("Should be $actDeg but is $corDeg")
        @req false "wrong b factor"
      end
    end
    res = res * bFactor
  end

  if check_degrees
    r = valency(R.gkm)
    ce = chern_number(e, R.gkm)
    corDeg = -(r-1) - d * ce
    if corDeg != _get_degree(res)
      println("Wrong h class for e=$e, d=$d")
      println("Deg is $(_get_degree(res)) but should be $corDeg")
      @req false "wrong h class"
    end
  end

  return res
end

# Calculate b(u,w,a) as in [Liu--Sheshmani, Lemma 4.5, p.16]. C is the coefficient ring
function _b(u::T, w::T, a::ZZRingElem) where T<:RingElem
  res = one(u) // one(u) # make sure this has FracFieldElem type (or QQFieldElem)
  if a >= 0
    for j in 0:a
      res = res // (w - j*u)
    end
  else 
    for j in 1:(-a-1)
      res = res * (w + j*u)
    end
  end
  return res
end

function GWTreeContribution(
  dt::GW_decorated_tree,
  P_input;
  check::Bool=true)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  # I assume this function is not used in fast mode, so we use equivariant parameters.
  C = dt.gkm.equivariantCohomology
  t = gens(C.coeffRing)

  res = Euler_inv(dt, t, C.edgeWeightClasses, C.pointEulerClasses)
  R = dt.gkm.equivariantCohomology
  con = get_connection(dt.gkm)

  # multiply by h classes
  for e in edges(dt.tree)
    res *= _h(imageOf(e, dt), edgeMult(e, dt), con, R, t, C.edgeWeightClasses; check)
  end

  #multiply by input class
  res *= Base.invokelatest(P_input.func, dt)

  return res
end

# Return the total h factor from the fiber direction.
function _h_VB(V::GKM_vector_bundle, e::Edge, d::Int, con::GKM_connection, R::GKM_cohomology_ring; check::Bool=true, check_degrees::Bool=false)::AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}

  gkm = con.gkm

  # Fast mode is not yet supported or vector bundles, so we use equivariant parameters.
  # TODO: unify compact and non-compact GKM graphs and make fast-mode available for both.
  C = R.coeffRing
  t = gens(C)

  if check
    @req con.gkm == R.gkm "GKM connection and cohomology ring don't belong to the same GKM graph"
    @req has_edge(gkm.g, e) "edge not found in GKM graph"
    @req d>0 "d is non-positive"
  end

  we = weight_class(e, R) # weight of the edge e

  # Start with the h factor from the base space.
  res = _h(e, d, con, R, t, R.edgeWeightClasses; check=check, check_degrees=check_degrees)

  # Apply h factors in fiber direction
  for i in 1:rank(V)

    wei = _fiber_summand_weight(src(e), i, V)
    ai = _fiber_connection_a(e, i, V)
    bFactor = _b(1//d * we, wei, d*ai)
    res = res * bFactor
  end

  return res
end