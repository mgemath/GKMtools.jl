@testset "Polynomial GKM spline ring" begin
  # This graph is used only through its GKM data; no module basis or ring
  # presentation is required.
  G = projective_space(GKMGraph, 1)
  H = polynomial_gkm_ring(G)
  S = equivariant_coefficient_ring(G)
  t = gens(S)
  e = first(edges(G))
  w = GKMtools.weight(G, e)
  alpha = sum(i -> w[i] * t[i], eachindex(t); init=zero(S))

  c = polynomial_class(G, [alpha, zero(S)])
  @test c isa GKMClass
  @test is_gkm_spline(c)
  @test restrictions(c) == [alpha, zero(S)]
  @test c * c == alpha * c
  @test c + zero(c) == c
  @test c^0 == one(c)
  @test localize_at_vertex(G, c, 1) == alpha

  localized = localize(c)
  @test is_gkm_spline(localized)
  @test is_gkm_spline(point_class(1, G))
  @test delocalize(G, localized) == c
  @test integrate(c) == integrate(localized)

  @testset "Scalar arithmetic" begin
    for class in (c, localized), scalar in (2, 2//3, 2.0, ZZ(2), QQ(2, 3))
      q = scalar isa AbstractFloat ? QQ(Rational{BigInt}(scalar)) : QQ(scalar)
      scaled = class * scalar
      divided = class / scalar
      @test scaled == scalar * class
      @test restrictions(scaled) == [x * q for x in restrictions(class)]
      @test restrictions(divided) == [x * inv(q) for x in restrictions(class)]
      @test class // scalar == divided
      @test divided * scalar == class
      @test parent(divided) === parent(class)
      @test graph(divided) === graph(class)
      @test typeof(divided) === typeof(class)
    end
    @test iszero(c * 0)
    @test_throws DivideError c / 0
    @test_throws DivideError c // 0
    @test_throws DivideError localized / 0
    @test_throws DivideError localized // 0
  end

  @test_throws ArgumentError polynomial_class(G, [one(S), zero(S)])

  rational = inv(gens_coeffRing(G)[1]) *
             one(localized)
  @test !is_gkm_spline(rational)
  @test_throws ArgumentError delocalize(G, rational)

  mktemp() do filename, io
    close(io)
    @test serialize_polynomial_class(filename, c) == filename
    restored = deserialize_polynomial_class(filename, G)
    @test restored == c
    @test restrictions(restored) == restrictions(c)
  end

  mktemp() do filename, io
    close(io)
    @test serialize_polynomial_class(filename, localized) == filename
    @test deserialize_polynomial_class(filename, G) == c
  end

  mktemp() do filename, io
    close(io)
    @test_throws ArgumentError serialize_polynomial_class(filename, rational)
  end
end
