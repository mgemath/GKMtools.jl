# Is this file deprecated? It is not called anywhere, hence I couldn't finish the update to flags. If we never use it, that is not a problem.

export Psi_pos_gen

function Psi_pos_gen(a, h)::EquivariantClass

  # rule = :(_Psi_pos_gen(dt, $a, $h))
  rule = :(_class_one(dt))
  return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, a, Int64[])
end


function _Psi_pos_gen(dt::Union{GW_decorated_tree, GW_decorated_graph}, a_input::Vector{Int64}, H::Dict{HodgeKey, QQFieldElem})

  ans = one(dt.gkm.equivariantCohomology.coeffRing)

  # The following line would only be needed for conversion to flags.
  # valG = valency(dt.gkm)

  for v in 1:n_vertices(dt.gkm.g)

    a = [a_input[i] for i in 1:length(dt.marks) if dt.marks[i] == v]  # exponent of the psi classes at vertex v

    # if isempty(a) || all(i -> i == 0, a)
    #   continue
    # end
  
    g = dt isa GW_decorated_tree ? 0 : dt.genus[v]; #= genus at vertex v=# println("g = $g")
    E_v = length(all_neighbors(dt.gkm.g, v)) # valency of vertex v
    Sum_ai = sum(a) # sum of exponents of psi classes at vertex v
    S_v = length(a) # number of markings at vertex v
  
    # check dimension constraints
    dim = 3*g - 3 + E_v + S_v
    (g > 0) && (Sum_ai > dim) && return zero(dt.gkm.equivariantCohomology.coeffRing)

    imV = imageOf(v, dt) # color of vertex v
    u = vcat((edgeMult(Edge(v, n), dt) .// weight_class(Edge(imV, imageOf(n, dt)), dt.gkm) for n in all_neighbors(dt.gkm.g, v))...,) #TODO: need to fix and adapt to flags.
    Sum_u = sum(u)
  
    if g == 0 # using the closed formula for genus 0
    
      (dim > -1) && (Sum_ai > dim) && return zero(dt.gkm.equivariantCohomology.coeffRing)

    #   temp = zero(dt.gkm.equivariantCohomology.coeffRing) # initialize answer
    
      if dim == -1 # necessary |S_v| == 1
        M = (-1)^a[1]
      else # E_v + n>2 and Sum_ai <= dim
        M = multinomial(dim - Sum_ai, a...,)
      end

    #   for n in all_neighbors(dt.g, v)
    #     temp += edgeMult(Edge(v, n), dt) // weight_class(Edge(imV, imageOf(n, dt)), dt.gkm)
    #   end
      
      ans *= M * (Sum_u^(-Sum_ai))
      continue
    end
  
    E_sigma_v = length(all_neighbors(dt.gkm.g, imV)) # valency of color of vertex v
    # u = [edgeMult(Edge(v, n), dt) // weight_class(Edge(imV, imageOf(n, dt)), dt.gkm) for n in all_neighbors(dt.g, v)]
    
    # TODO: The below is not yet flag-compatible. To implement it, one needs to add edge_weight_dict and t as arguments, as in the commented line below.
    w = [1 // weight_class(Edge(imV, n), dt.gkm) for n in all_neighbors(dt.gkm.g, imV)]
    # w = [1 // _flag_weight_class(dg.gkm, imV, i, t, edge_weight_dict) for i in 1:valG]

    

    NUMERATOR = zero(dt.gkm.equivariantCohomology.coeffRing)
    # DENOMINATOR = zero(dt.gkm.equivariantCohomology.coeffRing)

    for exponent in Oscar.weak_compositions(dim-Sum_ai, E_sigma_v + E_v)
    
      i = exponent[1:E_sigma_v] # exponent of the hodge classes
      any(i_a -> i_a > g, i) && continue # skip if any exponent of hodge class is greater than g
      lambda = [count(c -> c==t, i) for t in 1:g]

      j = exponent[(E_sigma_v+1):(E_sigma_v+E_v)] # exponent of the edge contributions

      COMMON = (-1)^(sum(i)) * prod(l -> w[l]^i[l], 1:E_sigma_v) * prod(k -> u[k]^j[k], 1:E_v)
      NUMERATOR   += hodge_integral(g, E_v + S_v, vcat(j, a), lambda, H) * COMMON
      println("NUMERATOR = $NUMERATOR")
    #   println("$g, $(E_v+S_v), $(vcat(j, a)), $lambda")
    end

    # for exponent in Oscar.weak_compositions(dim, E_sigma_v + E_v)
    
    #   i = exponent[1:E_sigma_v] # exponent of the hodge classes
    #   any(i_a -> i_a > g, i) && continue # skip if any exponent of hodge class is greater than g
    #   lambda = [count(c -> c==t, i) for t in 1:g]

    #   j = exponent[(E_sigma_v+1):(E_sigma_v+E_v)] # exponent of the edge contributions

    #   COMMON = (-1)^(sum(i)) * prod(l -> w[l]^i[l], 1:E_sigma_v) * prod(k -> u[k]^j[k], 1:E_v)
    #   # DENOMINATOR += hodge_integral(g, E_v, j, lambda, H) * (Sum_u^S_v) * COMMON
    #   DENOMINATOR += hodge_integral(g, E_v + S_v, vcat(j, zero(a)), lambda, H) * COMMON
    #   println("DENOMINATOR = $DENOMINATOR")
    # #   println("$g, $(E_v+S_v), $(vcat(j, zero(a))), $lambda")
    # end

    # iszero(DENOMINATOR) && return DENOMINATOR//one(dt.gkm.equivariantCohomology.coeffRing)
    # ans *= NUMERATOR//DENOMINATOR

    ans *= NUMERATOR
  end
println("ans in Psi= $ans")
  return ans
end