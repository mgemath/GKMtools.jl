using Test, Oscar, GKMtools

function qs_test_multiply(Q, a, b)
  result = Dict{Tuple{Int,Tuple{Vararg{Int}}},QQFieldElem}()
  for ((u,du),cu) in a, ((v,dv),cv) in b
    for ((w,dw),cw) in quantum_schubert_product(Q,u,v)
      d = Tuple(du[i]+dv[i]+dw[i] for i in eachindex(dw))
      key = (w,d)
      result[key] = get(result,key,QQ(0))+cu*cv*cw
    end
  end
  filter!(p -> !iszero(last(p)), result)
  return result
end

@testset "Quantum Schubert projective spaces" begin
  for n in 1:3
    Q = quantum_schubert_context(root_system(:A,n), collect(2:n))
    @test Q.novikov_degrees == [n+1]
    for u in 1:n+1, v in 1:n+1
      degree, codim = divrem(u+v-2,n+1)
      @test quantum_schubert_product(Q,u,v) == Dict((codim+1,(degree,)) => QQ(1))
    end
    @test quantum_schubert_product(Q,n+1,n+1;max_degree=0) == Dict()
  end
end

@testset "Quantum Schubert reconstruction and associativity" begin
  for (family, rank, levi) in [(:A,3,[1,3]), (:A,2,Int[]), (:B,2,[1]), (:B,2,[2]), (:G,2,[1]), (:G,2,[2])]
    Q = quantum_schubert_context(root_system(family,rank),levi)
    n = length(Q.representatives)
    table = quantum_schubert_table(Q)
    @test length(table) == div(n*(n+1),2)
    zd = Tuple(zeros(Int,length(Q.omitted_roots)))
    basis = [Dict((u,zd)=>QQ(1)) for u in 1:n]
    for i in Q.omitted_roots
      u = findfirst(==(reflection(simple_root(root_system(parent(Q.representatives[1])),i))), Q.representatives)
      for v in 1:n
        @test quantum_schubert_product(Q,u,v) == quantum_chevalley_product(Q,i,v)
      end
    end
    for u in 1:n, v in 1:n, w in 1:n
      @test qs_test_multiply(Q,quantum_schubert_product(Q,u,v),basis[w]) ==
        qs_test_multiply(Q,basis[u],quantum_schubert_product(Q,v,w))
    end
  end
end

@testset "Quantum Schubert graph conventions and localization" begin
  for (family, rank, levi) in [(:A,2,[2]), (:A,3,[1,3]), (:B,2,[1])]
    G = generalized_gkm_flag(root_system(family,rank),levi)
    Q = quantum_schubert_context(G)
    n = length(Q.representatives)
    order = sortperm(Q.codimensions)
    basis = [schubert_basis(G,v;representation=:polynomial) for v in order]
    S = GKMtools.equivariant_coefficient_ring(G)
    origin = fill(QQ(0),ngens(S))
    for i in 1:n, j in i:n
      cs = fill(zero(S),n)
      product = basis[i]*basis[j]
      for k in 1:n
        residual = product[order[k]]-sum((cs[l]*basis[l][order[k]] for l in 1:k-1);init=zero(S))
        ok, cs[k] = divides(residual,basis[k][order[k]])
        @test ok
        @test quantum_schubert_coefficient(Q,order[i],order[j],order[k]) == evaluate(cs[k],origin)
      end
    end
  end
  G = generalized_gkm_flag(root_system(:A,1),Int[])
  Q = quantum_schubert_context(G)
  beta = curve_class(G,Oscar.Edge(1,2))
  # Billey restrictions use the opposite root sign to the tangent weights.
  h = -schubert_basis(G,2)
  @test quantum_schubert_coefficient(Q,2,2,1,1) ==
    gromov_witten_nomarks(G,beta,[h,h,h];show_bar=false)
  G = generalized_gkm_flag(root_system(:A,2),[2])
  Q = quantum_schubert_context(G)
  beta = curve_class(G,Oscar.Edge(1,2))
  h,pt = -schubert_basis(G,2),schubert_basis(G,3)
  @test quantum_schubert_coefficient(Q,2,3,1,1) ==
    gromov_witten_nomarks(G,beta,[h,pt,pt];show_bar=false)
  G = generalized_gkm_flag(root_system(:A,3),[1,3])
  Q = quantum_schubert_context(G)
  a,b = findall(==(2),Q.codimensions)
  point = only(findall(==(4),Q.codimensions))
  beta = curve_class(G,Oscar.Edge(1,2))
  invariant = gromov_witten_nomarks(G,beta,
    [schubert_basis(G,a),schubert_basis(G,b),schubert_basis(G,point)];show_bar=false)
  @test invariant == 1
  @test quantum_schubert_coefficient(Q,a,b,Q.unit,1) == invariant
end

@testset "Quantum Schubert input validation and point" begin
  Q = quantum_schubert_context(root_system(:A,2),[2])
  @test_throws ArgumentError quantum_schubert_coefficient(Q,1,1,1,-1)
  @test_throws ArgumentError quantum_schubert_coefficient(Q,1,1,1,(1,2))
  @test_throws BoundsError quantum_schubert_product(Q,0,1)
  @test_throws ArgumentError quantum_chevalley_product(Q,2,1)
  @test_throws ArgumentError quantum_schubert_context(root_system(:A,2),[3])
  @test quantum_schubert_product(Q,"id","s1") == quantum_schubert_product(Q,1,2)
  P = quantum_schubert_context(root_system(:A,2),[1,2])
  @test quantum_schubert_product(P,1,1) == Dict((1,())=>QQ(1))
end

@testset "E7/P2 root-only reconstruction" begin
  Q = quantum_schubert_context(root_system(:E,7),[1,3,4,5,6,7])
  @test length(Q.representatives) == 576
  @test maximum(Q.codimensions) == 42
  @test Q.novikov_degrees == [14]
  divisor = only(findall(==(1),Q.codimensions))
  point = only(findall(==(42),Q.codimensions))
  # Compare a non-divisor product against D*(D*point), using only
  # Chevalley edges on the right-hand side (including q and q^2 terms).
  dsquared = quantum_chevalley_product(Q,2,divisor)
  @test length(dsquared) == 1
  (a,degree),scale = only(dsquared)
  @test degree == (0,)
  expected = Dict{Tuple{Int,Tuple{Vararg{Int}}},QQFieldElem}()
  for ((v,d),c) in quantum_chevalley_product(Q,2,point)
    for ((w,e),b) in quantum_chevalley_product(Q,2,v)
      key = (w,(d[1]+e[1],))
      expected[key] = get(expected,key,QQ(0))+c*b
    end
  end
  actual = Dict(k => scale*c for (k,c) in quantum_schubert_product(Q,a,point))
  @test actual == expected
  @test any(k -> k[2] == (2,),keys(actual))
end

@testset "Fixed-factor quantum Schubert serialization" begin
  Q = quantum_schubert_context(root_system(:A,2),[2])
  mktempdir() do dir
    path = joinpath(dir,"products.jls")
    saved = serialize_quantum_schubert_products(path,Q,3)
    loaded = GKMtools.Serialization.deserialize(path)
    @test loaded == saved
    @test loaded == [
      Dict((3,(0,)) => BigInt(1)),
      Dict((1,(1,)) => BigInt(1)),
      Dict((2,(1,)) => BigInt(1)),
    ]
    @test all(c isa BigInt for p in loaded for c in values(p))
    @test serialize_quantum_schubert_products(path,Q,"id") ==
      [Dict((v,(0,)) => BigInt(1)) for v in 1:3]
    @test_throws BoundsError serialize_quantum_schubert_products(path,Q,0)
    @test GKMtools.Serialization.deserialize(path) ==
      [Dict((v,(0,)) => BigInt(1)) for v in 1:3]
  end
end

@testset "Quantum Schubert multiplication matrices" begin
  Q = quantum_schubert_context(root_system(:A,2),[2])
  mktempdir() do dir
    products = serialize_quantum_schubert_products(joinpath(dir,"row.jls"),Q,2)
    M = quantum_schubert_matrix(products)
    S = base_ring(M)
    q = only(gens(S))
    @test M == matrix(S,[0 0 q; 1 0 0; 0 1 0])
    @test M^3 == q*identity_matrix(S,3)
    @test quantum_schubert_matrix(Q,2) == M
    @test quantum_schubert_matrix(products;at_q1=true) == matrix(QQ,[0 0 1; 1 0 0; 0 1 0])
    @test quantum_schubert_matrix(Q,2;at_q1=true) == quantum_schubert_matrix(products;at_q1=true)
  end
  products = [Dict((2,(1,0))=>BigInt(2),(2,(0,1))=>BigInt(3)),Dict((1,(0,0))=>BigInt(1))]
  M = quantum_schubert_matrix(products)
  q1,q2 = gens(base_ring(M))
  @test M[2,1] == 2q1+3q2
  @test M[1,2] == 1
  @test quantum_schubert_matrix(products;at_q1=true) == matrix(QQ,[0 1;5 0])
  @test quantum_schubert_matrix([Dict(),Dict()]) == zero_matrix(QQ,2,2)
  @test quantum_schubert_matrix([Dict((1,())=>BigInt(1))]) == identity_matrix(QQ,1)
  @test_throws ArgumentError quantum_schubert_matrix([Dict((2,(0,))=>1)])
  @test_throws ArgumentError quantum_schubert_matrix([Dict((1,(-1,))=>1)])
  @test_throws ArgumentError quantum_schubert_matrix([Dict((1,(0,))=>1),Dict((2,(0,0))=>1)])
end

@testset "Modular quantum Schubert arithmetic and bounded caches" begin
  for (family,r,levi) in [(:A,3,[1,3]),(:A,2,Int[]),(:B,2,[1]),(:G,2,[2])]
    R = root_system(family,r)
    Q = quantum_schubert_context(R,levi)
    exact = quantum_schubert_table(Q)
    for p in (11,101,1009)
      Qp = quantum_schubert_context(R,levi;p,max_cache_entries=40)
      F = Qp.field
      modular = quantum_schubert_table(Qp)
      for (uv,product) in exact
        reduced = Dict(k => F(c) for (k,c) in product if !iszero(F(c)))
        @test modular[uv] == reduced
      end
      @test length(Qp.cache) <= 40
      @test length(Qp.diagonal_cache) <= 40
      @test length(Qp.restriction_cache) <= 40
      @test length(Set(Qp.diagonal)) == length(Qp.representatives)
      @test all(!iszero,Qp.self_restrictions)
      @test Qp.prime == p
      for u in eachindex(Qp.representatives), w in eachindex(Qp.representatives)
        @test GKMtools._qs_below(Qp,u,w) ==
          (u == w || Qp.representatives[u] < Qp.representatives[w])
      end
    end
  end
  R = root_system(:A,3)
  unlimited = quantum_schubert_context(R,[1,3])
  tiny = quantum_schubert_context(R,[1,3];max_cache_entries=5)
  @test quantum_schubert_table(tiny) == quantum_schubert_table(unlimited)
  @test length(tiny.cache) <= 5
  @test length(tiny.diagonal_cache) <= 5
  @test length(tiny.restriction_cache) <= 5
  @test_throws ArgumentError quantum_schubert_context(R,[1,3];p=15)
  @test_throws ArgumentError quantum_schubert_context(R,[1,3];p=5)
  @test_throws ArgumentError quantum_schubert_context(R,[1,3];p=101.0)
  @test_throws ArgumentError quantum_schubert_context(R,[1,3];max_cache_entries=-1)
  Qp = quantum_schubert_context(root_system(:A,2),[2];p=101)
  mktempdir() do dir
    path = joinpath(dir,"modular.jls")
    products = serialize_quantum_schubert_products(path,Qp,3)
    @test products == GKMtools.Serialization.deserialize(path)
    @test all(c isa BigInt && 0 <= c < 101 for row in products for c in values(row))
    @test quantum_schubert_matrix(Qp,3;at_q1=true) ==
      quantum_schubert_matrix(products;p=101,at_q1=true)
    @test characteristic(base_ring(quantum_schubert_matrix(products;p=101,at_q1=true))) == 101
    @test quantum_schubert_matrix(Qp,3) == quantum_schubert_matrix(products;p=101)
  end
  G = generalized_gkm_flag(root_system(:A,2),[2])
  @test quantum_schubert_table(quantum_schubert_context(G;p=101)) == quantum_schubert_table(Qp)
end

@testset "Large prime and E7/P5 modular reconstruction" begin
  p = BigInt(next_prime(ZZ(typemax(Int))))
  Q = quantum_schubert_context(root_system(:A,2),[2];p)
  @test GKMtools._qs_products(Q,3) == [
    Dict((3,(0,))=>BigInt(1)), Dict((1,(1,))=>BigInt(1)), Dict((2,(1,))=>BigInt(1))]
  Q = quantum_schubert_context(root_system(:E,7),[1,2,3,4,6,7];p=1_000_000_007)
  @test length(Q.representatives) == 4032
  @test Q.novikov_degrees == [10]
  @test maximum(Q.codimensions) == 50
  divisor = only(Q.degree_vertices[2])
  point = only(Q.degree_vertices[end])
  dsquared = quantum_chevalley_product(Q,5,divisor)
  expected = Dict{Tuple{Int,Tuple{Vararg{Int}}},elem_type(Q.field)}()
  for ((v,d),c) in quantum_chevalley_product(Q,5,point)
    for ((w,e),b) in quantum_chevalley_product(Q,5,v)
      key = (w,(d[1]+e[1],))
      expected[key] = get(expected,key,zero(Q.field))+c*b
    end
  end
  actual = empty(expected)
  for ((u,d),c) in dsquared
    @test d == (0,)
    for (k,b) in quantum_schubert_product(Q,u,point)
      actual[k] = get(actual,k,zero(Q.field))+c*b
    end
  end
  @test actual == expected
  @test length(Q.cache) <= Q.max_cache_entries
end
