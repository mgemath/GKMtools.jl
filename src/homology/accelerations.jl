# The compact Morse-index formula is used under the compact Hamiltonian
# hypothesis of HM26 §4.2. It must never be applied to standalone flags.
_homology_is_compact(data::GKMCombinatorialData) = all(
    v -> length(flags(data, v)) == length(all_neighbors(graph(data), v)),
    1:nv(graph(data)),
)

# Keep orbifold normalization on the full relation algorithm.
_homology_use_compact_formula(data::GKMCombinatorialData) = _homology_is_compact(data)
_homology_use_compact_formula(::GKMCombinatorialData{R,V,F}) where {R,V,F<:AbstractOrbifoldFlagWeight} = false

# Surjective integral projection: the first two columns have determinant one.
# Moment-curve searches use ZZ even for their parameters, avoiding overflow.
function _homology_rank_two_projection(data::GKMCombinatorialData)
    r = rank_torus(data)
    r <= 2 && return identity_matrix(QQ, r)
    ws = [[QQFieldElem[weight(data, v, i)[j] for j in 1:r]
           for i in eachindex(flags(data, v))] for v in 1:nv(graph(data))]
    @req all(w -> any(!iszero, w), Iterators.flatten(ws)) "Flag weights must be nonzero"
    @req all(length(weights) < 2 || all(
        rank(matrix(QQ, [weights[i], weights[j]])) == 2
        for i in eachindex(weights) for j in (i+1):length(weights)
    ) for weights in ws) "Flag weights must be pairwise independent"
    t = ZZ(0)
    a = QQFieldElem[]
    av = Vector{QQFieldElem}[]
    while true
        a = QQFieldElem[t^(j-1) for j in 1:r]
        av = [[dot(a,w) for w in weights] for weights in ws]
        all(!iszero, Iterators.flatten(av)) && break
        t += 1
    end
    s = ZZ(0)
    while true
        b = QQFieldElem[j == 1 ? 0 : s^(j-2) for j in 1:r]
        bv = [[dot(b,w) for w in weights] for weights in ws]
        if all(all(av[v][i]*bv[v][j] != av[v][j]*bv[v][i]
                   for i in eachindex(ws[v]) for j in (i+1):length(ws[v]))
               for v in eachindex(ws))
            return matrix(QQ, [a,b])
        end
        s += 1
    end
end

# Values of a nonzero functional on a known rank-one edge lattice determine
# its primitive quotient. They need not themselves be primitive or positive.
function _primitive_homology_quotient(edge_lattice, values)
    denominator_lcm = foldl(lcm, (denominator(QQ(x)) for x in values); init=ZZ(1))
    integers = ZZRingElem[ZZ(denominator_lcm * x) for x in values]
    divisor = foldl(gcd, integers; init=ZZ(0))
    @req !iszero(divisor) "The homology functional must be nonzero"
    first_nonzero = findfirst(!iszero, integers)
    divisor *= sign(integers[first_nonzero])
    H2 = free_module(ZZ, 1)
    quotient = ModuleHomomorphism(edge_lattice, H2,
        [divexact(x, divisor)*gen(H2,1) for x in integers])
    return H2, quotient
end

function _homology_relations(data, edge_list, edge_to_gen; project::Bool)
    cycles = _fundamental_cycles(data, edge_list, edge_to_gen)
    projection = project ? _homology_rank_two_projection(data) : identity_matrix(QQ, rank_torus(data))
    weights = Dict(e => projection * matrix(QQ, rank_torus(data), 1, _homology_weight(data,e))
                   for edge in edge_list for e in (edge, reverse(edge)))
    ncoords = nrows(projection)
    relations = zero_matrix(QQ, Base.checked_mul(ncoords, length(cycles)), length(edge_list))
    for (i, cycle) in enumerate(cycles), j in 1:ncoords, e in cycle
        relations[(i-1)*ncoords+j, edge_to_gen[e]] += weights[e][j,1]
    end
    denominator_lcm = foldl(lcm, (denominator(x) for x in relations); init=ZZ(1))
    return change_base_ring(ZZ, denominator_lcm * relations)
end
