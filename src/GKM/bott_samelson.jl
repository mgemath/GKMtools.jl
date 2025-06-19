export bott_samelson

@doc raw"""
    bott_samelson(S::Vector{RootSpaceElem}; printEdges::Bool = true)

We construct the Bott-Samelson desingularization of a Schubert variety given by `S`.
"""
function bott_samelson(S::Vector{RootSpaceElem}; printEdges::Bool = true)
  @req all(s -> is_simple_root(s), S) "Roots for Bott-Samelson variety must be simple."
  @req length(S) > 0 "Neet at least one root for Bott-Samelson variety."
  R = root_system(S[1])
  @req all(s -> root_system(s) == R, S) "Roots don't belong to the same root system."
  @req length(S) <= 16 "At most 16 weights are allowed for Bott-Samelson varieties currently."

  nr = length(S)
  refls = [reflection(s) for s in S]

  # these will be all partial products of refls in order.
  epsilons = Vector{WeylGroupElem}(undef, 2^(nr-1))
  epsilons[1] = one(parent(refls[1]))

  edgeFlags = BitVector(undef, 2^(2*nr))
  nEdges = 0
  #edgeFlags = BitVector(undef, 2^(2*nr))
  fatEdges = Vector{UInt32}()
  
  for n in 1:nr
    w = refls[n]
    n > 1 && (epsilons[(2^(n-2)+1):(2^(n-1))] = epsilons[1:2^(n-2)] .* w)
    #println("epsilons = $epsilons")

    rMax = UInt16(2^(n-1)-1) # bitstrings of length at most n-1
    for x in UInt16(0):rMax
      for y in x:rMax

        x0 = _bs_x0(x, n)
        x1 = _bs_x1(x, n)

        if x == y
          edgeFlags[_bs_edge_index(x0, x1, nr)] = true # add edge (0x, 1x): Note, we always use (x < y) for edges!
          nEdges += 1
          continue
        end

        y0 = _bs_x0(y, n)
        y1 = _bs_x1(y, n)
        x0y1 = _bs_edge_index(x0, y1, nr)
        x1y0 = _bs_edge_index(y0, x1, nr)
        x1y1 = _bs_edge_index(x1, y1, nr)

        if !edgeFlags[_bs_edge_index(x, y, nr)]
          edgeFlags[x0y1] = false
          edgeFlags[x1y1] = false
          edgeFlags[x1y0] = false
          continue
        end

        
        k = _get_first_differing_bit(x, y)
        t = _get_bs_type(x, y, n, k, S, epsilons)

        if t == 1
          edgeFlags[x1y1] = true # note x1 < y1 since x < y etc.
          edgeFlags[x0y1] = false
          edgeFlags[x1y0] = false
          nEdges += 1
        elseif t == 2
          edgeFlags[x1y1] = false
          edgeFlags[x0y1] = true
          edgeFlags[x1y0] = true
          nEdges += 2
        elseif t == 3
          edgeFlags[x1y1] = true
          edgeFlags[x0y1] = true
          edgeFlags[x1y0] = false
          nEdges += 2
        elseif t == 4
          edgeFlags[x1y1] = true
          edgeFlags[x0y1] = false
          edgeFlags[x1y0] = true
          nEdges += 2
        else
          @req false "_get_bs_type returned value outside 1:4. We got t=$t"
        end

        # in final round, record the fat edges corresponding to 1-dimensional families of T-stable curves.
        if n == nr
          if t == 2
            push!(fatEdges, _bs_edge_index(x0, y0, nr))
          elseif t == 3
            push!(fatEdges, _bs_edge_index(x0, y1, nr))
          elseif t == 4
            push!(fatEdges, _bs_edge_index(y0, x1, nr))
          end
        end
      end
    end
  end
  
  println("Bott samelson graph in type $R defined by $S")
  println("Is GKM: $(isempty(fatEdges))")
  println("Number of vertices: $(2^nr)")
  println("Number of edges: $nEdges")
  println("Number of fat edges: $(length(fatEdges))")
  
  !printEdges && return

  println("Edges:")

  nEdges = 0
  rMax = UInt16(2^nr-1)
  for x in UInt16(0):rMax
    xLabel = reverse(bitstring(x)[(16-nr+1):16])
    for y in UInt16(x+1):rMax

      edgeIndex = _bs_edge_index(x, y, nr)
      !edgeFlags[edgeIndex] && continue
      nEdges += 1
      yLabel = reverse(bitstring(y)[(16-nr+1):16])
      if edgeIndex in fatEdges
        println("$xLabel -> $yLabel (fat)")
      else
        println("$xLabel -> $yLabel")
      end
    end
  end
  println("Number of edges: $nEdges")

end

# This is always used with x < y
function _bs_edge_index(x::UInt16, y::UInt16, nr::Int64)::UInt32
  return (UInt32(x) << nr) + y
end

# Return 1, 2, 3, or 4, according to [Withrow, Thm 1]
function _get_bs_type(x::UInt16, y::UInt16, n::Int64, k::Int64, S::Vector{RootSpaceElem}, epsilons::Vector{WeylGroupElem})
  ck = UInt16(2^k-1) # use _ & ck to truncate to the first k bits.
  x0 = _bs_x0(x, n)
  y0 = _bs_x0(y, n)
  ex0k = S[k] * inv(epsilons[(x0 & ck) + 1])
  ex0n = S[n] * inv(epsilons[x0 + 1])
  ey0k = S[k] * inv(epsilons[(y0 & ck) + 1])
  ey0n = S[n] * inv(epsilons[y0 + 1])

  # println("x=$x, y=$y, n=$n, k=$k: x0=$x0, y0=$y0, ex0k=$ex0k, ex0n=$ex0n, ey0k=$ey0k, ey0n=$ey0n, ck=$ck")
  # println("ex0k = $(S[k]) * $(inv(epsilons[(x0 & ck) + 1]))")
  # println("ex0n = $(S[n]) * $(inv(epsilons[ x0 + 1]))")
  # println("ey0k = $(S[k]) * $(inv(epsilons[(y0 & ck) + 1]))")
  # println("ey0n = $(S[n]) * $(inv(epsilons[y0 + 1]))")

  (ex0k != ex0n) && (ey0k != ey0n) && return 1
  (ex0k == ex0n) && (ey0k == ey0n) && return 2
  (ex0k == ex0n) && (ey0k != ey0n) && return 3
  (ex0k != ex0n) && (ey0k == ey0n) && return 4
end

function _bs_x0(x::UInt16, n::Int64)::UInt16
  return x & UInt16(2^(n-1) - 1)
end

function _bs_x1(x::UInt16, n::Int64)::UInt16
  return _bs_x0(x, n) + UInt16(2^(n-1))
end

function _get_first_differing_bit(x::UInt16, y::UInt16)
  b = UInt16(1)
  for n in 0:15
    if b & (x >> n) != b & (y >> n)
      return n+1
    end
  end
  return nothing # this happens iff x == y
end