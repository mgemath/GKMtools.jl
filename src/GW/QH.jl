export QH_print_structure_constants_in_basis, QH_print_structure_constants

@doc raw"""
    QH_structure_constants(G::AbstractGKM_graph; refresh::Bool=false)

Return the structure constants of the equivariant quantum cohomology $QH_T^*(X)$ where $X$ is the GKM space realizing the GKM graph.
The structure constants are returned with respect to the [fixed point basis](#The-fixed-point-basis) as input basis and the [standard basis](#The-standard-basis) as output basis
(see [Mixing the bases](#Mixing-the-bases)).
This has the advantage that all entries are polynomials in the equivariant parameters when `G` is the GKM graph of a compact Hamiltonian or projective GKM space.

!!! warning
    - This requires `is_strictly_nef(G)==true`, as this guarantees that there are at most finitely many curve classes $\beta$ with non-zero coefficients for $q^\beta$.
    - If `is_strictly_nef(G)==false`, use the method of [`QH_structure_constants`](@ref) below that specifies a specific $\beta$.

!!! note
    - As this computation might be expensive, the result is stored in `G` for later use. If the requested structure constants have been computed before, they will not be
      computed afresh unless the optional argument `refresh` is set to `true`.

# Output format:
The output type is a Julia dictionary of type `Dict{CurveClass_type, Array{Any, 3}}`. For each curve class $\beta$, the corresponding value is a 3-dimensional array of type `Array{Any, 3}`.
If `ans` denotes this dictionary, then
`ans[beta][i, j, k]` is the the $q^\beta$-coefficient of $PD(v_i) \ast PD(v_j)$ localized at $v_k$, where $v_i, v_j, v_k$ represent the fixed points with indices $i,j,k$, respectively, 
and $PD$ represents the equivariant Poincaré dual (cf. [fixed point basis](#The-fixed-point-basis)).

# Optional arguments:
 - `refresh::Bool`: `false` by default. If `true`, then this will overwrite any previously calculated $QH_T$ structure constants of `G`.

# Example
In the following example, we calculate the structure constants of $\mathbb{P}^1$.
The fixed point basis $f_1,f_2$ (see above) is given by `f1 = (t1 - t2)*e[1]` and `f2 = (-t1 + t2)*e[2]`. 
Let us check that
```math
f_1 \ast f_1 = (t_1 - t_2) f_1 + q = (t_1 - t_2)^2 e_1 + q(e_1 + e_2)
```
noting that $e_1+e_2=1$.
Similarly, we will check that
```math
f_1 \ast f_2 = q = q(e_1 + e_2)
```
and
```math
f_2 \ast f_2 = (-t_1 + t_2) f_2 + q = (-t_1 + t_2)^2 e_2 + q(e_1 + e_2).
```
The following code verifies these identities.

```jldoctest QH_structure_constants_all
julia> P1 = projective_space(GKM_graph, 1);

julia> S = QH_structure_constants(P1; show_progress=false);

julia> curve_classes = collect(keys(S));

julia> S[curve_classes[1]] # curve class 0
2×2×2 Array{Any, 3}:
[:, :, 1] =
 t1^2 - 2*t1*t2 + t2^2  0
 0                      0

[:, :, 2] =
 0  0
 0  t1^2 - 2*t1*t2 + t2^2

julia> S[curve_classes[2]] # line class
2×2×2 Array{Any, 3}:
[:, :, 1] =
 1  1
 1  1

[:, :, 2] =
 1  1
 1  1
```
"""
function QH_structure_constants(G::AbstractGKM_graph; refresh::Bool=false, show_progress::Bool=true)

  @req is_strictly_nef(G) "G is not strictly NEF, so need to specify beta"

  if !refresh && G.know_all_QH_structure_consts
    return G.QH_structure_consts
  end

  if refresh
    G.QH_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()
  end

  # Max deg of cohomology class (primitive wrt equivariant parameters): 2n
  # Max deg of product of two such: 4n
  # Degree of q^beta: 2c1(beta)
  # Thus, assuming G is strictly NEF, we only need 0 <= c1(beta) <= 2n
  n = valency(G)
  maxChernNumber = 2*n

  # Generating P_input takes a little if we have many vertices. So we calculate it here and reuse it.
  show_progress && println("Building classes to integrate over the moduli space...")
  nv = n_vertices(G.g)
  evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
  P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
  show_progress && println("Starting integration:")

  for c in 0:maxChernNumber
    for beta in _effectiveClassesWithChernNumber(GKM_second_homology(G), c)
      show_progress && println("Calculating structure constants for $beta, Chern number $c:")
      QH_structure_constants(G, beta; refresh = refresh, P_input = P_input, show_progress = show_progress)
    end
  end

  G.know_all_QH_structure_consts = true

  return G.QH_structure_consts
end

@doc raw"""
    QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false, P_input=nothing, show_progress::Bool=true)

Return the $q^\beta$-coefficients of the structure constants of the equivariant quantum cohomology $QH_T^*(X)$, where $X$ is the GKM space (see [Definition](../GKM/GKM.md#Definition)) realizing the GKM graph.

!!! note
    - As this computation might be expensive, the result is stored in `G` for later use. If the requested structure constants have been computed before, they will not be
      computed afresh unless the optional argument `refresh` is set to `true`.

# Output format:
The output type is `Array{Any, 3}`.
If `ans` denotes the returned object, then
`ans[i, j, k]` is the the $q^\beta$-coefficient of $PD(v_i) \ast PD(v_j)$ localized at $v_k$, where $v_i, v_j, v_k$ represent the fixed points with indices $i,j,k$, respectively, 
and $PD$ represents the Poincaré dual.

# Optional arguments:
 - `refresh::Bool`: `false` by default. If `true`, then this will overwrite any previously calculated $QH_T$ structure constants of `G`.

# Example

Let us repeat the example of $\mathbb{P}^1$ above by specifying each curve class $\beta$ separately.

```jldoctest QH_structure_constants_edge
julia> P1 = projective_space(GKM_graph, 1);

julia> beta = curve_class(P1, Edge(1, 2));

julia> QH_structure_constants(P1, 0*beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 t1^2 - 2*t1*t2 + t2^2  0
 0                      0

[:, :, 2] =
 0  0
 0  t1^2 - 2*t1*t2 + t2^2

julia> QH_structure_constants(P1, beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 1  1
 1  1

[:, :, 2] =
 1  1
 1  1

julia> QH_structure_constants(P1, 2*beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 0  0
 0  0

[:, :, 2] =
 0  0
 0  0

julia> QH_structure_constants(P1, -1 * beta; show_progress=false)
2×2×2 Array{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}, 3}:
[:, :, 1] =
 0  0
 0  0

[:, :, 2] =
 0  0
 0  0
```
"""
function QH_structure_constants(G::AbstractGKM_graph, beta::CurveClass_type; refresh::Bool=false, P_input=nothing, show_progress::Bool=true)

  if refresh || !haskey(G.QH_structure_consts, beta)

    R = G.equivariantCohomology.coeffRingLocalized
    nv = n_vertices(G.g)
    res = zeros(R, nv, nv, nv)
  
    if beta == zero(parent(beta))
      for i in 1:nv
        res[i, i, i] = (euler_class(i, G)//1)^2
      end
      G.QH_structure_consts[beta] = res
      return res
    end
    if !is_effective(G, beta)
      G.QH_structure_consts[beta] = res
      return res
    end
  
    # If we calculate for multiple beta, it makes sense to reuse P_input, hence the optional argument.
    # Generating P_input actually takes a surprisingly long amount of time for growing number of vertices.
    # For example for n=16 vertices it takes around 10 seconds.
    if isnothing(P_input)
      show_progress && println("Building classes to integrate over the moduli space...")
      evClasses = [ev(i, point_class(j, G)) for i in 1:3, j in 1:nv]
      P_input = [evClasses[1, i]*evClasses[2, j]*evClasses[3, k] for i in 1:nv, j in 1:nv, k in 1:nv]
      show_progress && println("Starting integration:")
    end
    res = gromov_witten(G, beta, 3, P_input; show_bar=show_progress)
    G.QH_structure_consts[beta] = res
    return res
  end
  return G.QH_structure_consts[beta]
end

@doc raw"""
    quantum_product(G::AbstractGKM_graph, beta::CurveClass_type, class1, class2; useStructureConstants::Bool = true, fastMode::Bool = false, distantVertex::Int64 = 1, twist_class::Union{Nothing, EquivariantClass}=nothing, show_progress::Bool=true)

Calculate the $q^\beta$-coefficient of the equivariant quantum product of the equivariant cohomology classes `class1` and `class2` on `G`.

If the optional argument `useStructureConstants` is set to `false`, then this will always calculate the relevant Gromov--Witten invariants
freshly using [`gromov_witten`](@ref), even if they have been calculated before.

# Fast mode
The optional argument `fastMode` must only be set to `true` when one is certain that the output is a degree zero cohomology class, i.e. a rational number.
It is not yet supported in combination with `useStructureConstants=true`.
When fast mode is used, the result is calculated using a single 3-point Gromov--Witten invariant, which takes
as arguments `class1`, `class2`, and the point class of the vertex with number `distantVertex`.
By default, the optional argument `distantVertex` is `1`, and its value does not change the result (if the result is known to be a rational number).
When `class1` and `class2` are also point classes, performance may be optimized by picking a value for
`distantVertex` such that there are relatively few trees in curve class `beta` meeting the points given by `class1`, `class2`, and `distantVertex`.

# Twisting
The experimental optional argument `twist_class` allows to multiply by any additional `EquivariantClass` before integrating over $\overline{\mathcal{M}}_{0,3}(X;\beta)$ in
the definition of the quantum product.
Currently, it can only be used in combination with `useStructureConstants=false`.
See [`reduced_virtual_zero_section`](@ref) for a twisted example.

# Progress bar
The optional progress bars for each underlying call to [`gromov_witten`](@ref) can be activated by setting `show_progress=true`.

# Example

Let us compute the equivariant quantum products of some point classes on $\mathbb{P}^2$.

```jldoctest quantum_product
julia> P2 = projective_space(GKM_graph, 2);

julia> beta = curve_class(P2, Edge(1, 2));

julia> quantum_product(P2, beta, point_class(P2, 1), point_class(P2, 2))
(t1 - t3, t2 - t3, 0)

julia> quantum_product(P2, 0*beta, point_class(P2, 1), point_class(P2, 2))
0

julia> quantum_product(P2, 2*beta, point_class(P2, 1), point_class(P2, 2))
(0, 0, 0)
```
"""
function quantum_product(
  G::AbstractGKM_graph,
  beta::CurveClass_type,
  class1, #class1 and class2 should be free module elements over the coefficient ring (or its frac field)
  class2;
  useStructureConstants::Bool = true,
  fastMode::Bool = false,
  distantVertex::Int64 = 1,
  twist_class::Union{Nothing, EquivariantClass}=nothing,
  show_progress::Bool=false
)

  @req !(useStructureConstants && fastMode) "Fast mode and structure constants are not simultaneously supported yet."
  @req !(useStructureConstants && !isnothing(twist_class)) "Twisting and structure constants are not simultaneously supported yet." #TODO: implement structure constants for twisting.

  if iszero(beta)

    # Make sure class1 and class2 are elements of the free module with the same coefficient ring.
    # (Either coefficient ring could be a polynomial ring or its fraction field)
    class1_is_fraction = base_ring(class1) == fraction_field(base_ring(class1))
    class2_is_fraction = base_ring(class2) == fraction_field(base_ring(class2))

    if class1_is_fraction && !class2_is_fraction
      # need to convert class2 to fraction field coefficients
      R = parent(class1)
      class2 = R([class2[i] for i in 1:rank(R)])
    elseif class2_is_fraction && !class1_is_fraction
      # need to convert class1 to fraction field coefficients
      R = parent(class2)
      class1 = R([class1[i] for i in 1:rank(R)])
    end

    return class1 * class2
  end

  nv = n_vertices(G.g)

  if fastMode
    @req distantVertex > 0 && distantVertex <= nv "distantVertex must be in 1:nv"
    GW_invt = gromov_witten(G, beta, 3, ev(1, class1) * ev(2, class2) * ev(3, point_class(distantVertex, G)); fast_mode=true, show_bar=show_progress)
    return GW_invt * one(G.equivariantCohomology)
  end

  if useStructureConstants
    C = QH_structure_constants(G, beta; show_progress=show_progress)
    res = zero(G.equivariantCohomology.cohomRingLocalized)
    for i in 1:nv, j in 1:nv, k in 1:nv
      eulerI = euler_class(i, G)
      eulerJ = euler_class(j, G)
      res += (C[i, j, k] * class1[i] * class2[j] // eulerI // eulerJ) * gens(G.equivariantCohomology.cohomRingLocalized)[k]
    end
    return res
  end

  P_input = [ev(1, class1)*ev(2, class2)*ev(3, point_class(v, G)) for v in 1:nv]
  if !isnothing(twist_class)
    P_input = [twist_class * p for p in P_input] # twist_class .* P_input errors for some reason...
  end
  GW_invts = gromov_witten(G, beta, 3, P_input; show_bar=show_progress)

  return sum([GW_invts[i] * gens(G.equivariantCohomology.cohomRingLocalized)[i] for i in 1:nv])
end


@doc raw"""
    QH_structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix; setPreferredBasis::Bool=false)

Return all structure constants of `G` that have been calculated so far with respect to the given basis.
A smart choice of basis can drastically simplify the presentation of the ring $QH_T^*(X)$.
Here, we use the given basis `b` both as input basis and as output basis for computing the structure constants (cf. [mixing the bases](#Mixing-the-bases)).

!!! note
    - This does not calculate any structure constants afresh. To do so, use [`QH_structure_constants`](@ref).
    - This will omit any curve classes in which all structure constants are zero.

# Output format:
The same as that of [`QH_structure_constants`](@ref), i.e. of type `Dict{CurveClass_type, Array{Any, 3}}`.

# Arguments
 - `G::AbstractGKM_graph`: The GKM graph whose quantum cohomology is of interest.
 - `b::Matrix`: A matrix whose rows are the desired $H_T^*(\text{pt};\mathbb{Q})$-linear basis of $H_T^*(X;\mathbb{Q})$.
    The element `b[i,j]` is the localization to the `j`-th fixed point of the `i`-th basis element.
 - `setPreferredBasis::Bool=false`: Optional argument. If set to `true`, all future quantum cohomology classes of this space
    will be printed with respect to the given base.

# Examples
This example shows that $QH_T(\mathbb{P}^1;\mathbb{Q}) \cong\mathbb{Q}[t_1, t_2, e,q]/(e^2 - (t_1-t_2)e - q)$ where $e=PD([1:0])$ and $q$
corresponds to the curve class $[\mathbb{P}^1]\in H_2(\mathbb{P}^1;\mathbb{Z})$.
```jldoctest QH_structure_constants_in_basis
julia> P1 = projective_space(GKM_graph, 1);

julia> QH_structure_constants(P1; show_progress=false);

julia> P1 = projective_space(GKM_graph, 1);

julia> QH_structure_constants(P1; show_progress=false)
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Array{Any, 3}} with 2 entries:
  (0) => [t1^2 - 2*t1*t2 + t2^2 0; 0 0;;; 0 0; 0 t1^2 - 2*t1*t2 + t2^2]
  (1) => [1 1; 1 1;;; 1 1; 1 1]

julia> t1, t2 = gens(P1.equivariantCohomology.coeffRing);

julia> base = [1 1; t1-t2 0];

julia> QH_structure_constants_in_basis(P1, base)
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Array{Any, 3}} with 2 entries:
  (0) => [1 0; 0 0;;; 0 1; 1 t1 - t2]
  (1) => [0 0; 0 1;;; 0 0; 0 0]
```

Similarly, choosing a nice basis simplifies the presentation of $QH_T(\mathbb{P}^2)$.
By the below, it is isomorphic as $H_T^*(\text{pt};\mathbb{Q})$-algebra to
$\mathbb{Q}[t_1,t_2,t_3, e, q]/(e(e-t_1+t_2)(e-t_1+t_3) - q)$, 
where $e = PD(\mathbb{P}^1_{[x:y:0]})$.
```jldoctest QH_structure_constants_in_basis
julia> P2 = projective_space(GKM_graph, 2);

julia> QH_structure_constants(P2; show_progress=false); # Calculate all relevant structure constants.

julia> (t1, t2, t3) = gens(P2.equivariantCohomology.coeffRing);

julia> base = [1 1 1; t1-t3 t2-t3 0 ; (t1-t2)*(t1-t3) 0 0];

julia> S = QH_structure_constants_in_basis(P2, base);

julia> beta = curve_class(P2, Edge(1, 2));

julia> S[beta][:,:,1]
3×3 Matrix{Any}:
 0  0  0
 0  0  1
 0  1  t1 - t2

julia> S[beta][:,:,2]
3×3 Matrix{Any}:
 0  0  0
 0  0  0
 0  0  1

julia> S[beta][:,:,3]
3×3 Matrix{Any}:
 0  0  0
 0  0  0
 0  0  0
```

"""
function QH_structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix; setPreferredBasis::Bool=false)
  s = size(b)
  @req s[1] == s[2] "Base matrix must be square"
  @req s[1] == n_vertices(G.g) "dimension of basis elements must be number of vertices of GKM graph"

  FF = G.equivariantCohomology.coeffRingLocalized
  nv = n_vertices(G.g)
  bMatrix = zero_matrix(FF, nv, nv)
  for k in keys(b)
    bMatrix[k] = FF(b[k])
  end
  @req is_invertible(bMatrix) "Base matrix must be invertible over fraction field"
  bMatrixInv = inv(transpose(bMatrix))
  if setPreferredBasis
    G.QH_preferred_basis = bMatrixInv
  end

  res = Dict{CurveClass_type, Array{Any, 3}}()
  g = gens(G.equivariantCohomology.cohomRingLocalized)
  baseClasses = [sum([bMatrix[i,j] * g[j] for j in 1:nv]) for i in 1:nv]
  for beta in keys(G.QH_structure_consts)
    resBeta = zeros(G.equivariantCohomology.coeffRingLocalized, (nv, nv, nv)...)
    for i in 1:nv, j in 1:nv
      prodIJ = quantum_product(G, beta, baseClasses[i], baseClasses[j])
      prodIJVect = [prodIJ[k] for k in 1:nv]
      inBasis = bMatrixInv * prodIJVect
      resBeta[i,j,:] = inBasis
    end
    all(x -> iszero(x), resBeta) && continue
    res[beta] = resBeta
  end
  return res
end

function QH_print_structure_constants(G::AbstractGKM_graph)
  nv = n_vertices(G.g)
  l = ["e_{$i}" for i in 1:nv]
  b = Matrix{Any}(undef, nv, nv)
  for i in 1:nv
    b[i, :] = repeat([0], nv)
    b[i, i] = euler_class(i, G)
  end
  QH_print_structure_constants_in_basis(G, b, l)
end

function QH_print_structure_constants_in_basis(G::AbstractGKM_graph, b::Matrix, l::Vector{String})
  S = QH_structure_constants_in_basis(G, b)
  betas = keys(S)
  nv = n_vertices(G.g)
  for i in 1:nv
    for j in i:nv
      printedIJ = false
      print(l[i] * " * " * l[j] * " &=& ")
      for b in betas
        betaTerms = count(x -> !iszero(x), S[b][i, j, :])
        printedIJ && betaTerms > 0 && print(" + ")
        betaTerms > 1 && print("(")
        printedBeta = false
        for k in 1:nv
          p = S[b][i, j, k]
          iszero(p) && continue
          printedBeta && print(" + ")
          _print_coefficient(p)
          print(" " * l[k])
          printedBeta = true
          printedIJ = true
        end
        if printedBeta
          betaTerms > 1 && print(")")
          _print_curve_q(b)
        end
      end
      !printedIJ && print("0")
      print("\\\\\n")
    end
  end
end

function _print_curve_q(b)
  #TODO: add custom labels for curve classes
  if !iszero(b)
    print(" q^{$(b)}")
  end
end

function _print_coefficient(s)
  if iszero(s)
    print("0")
    return
  elseif s == 1
    return
  else
    if denominator(s) == 1
      _factor_and_print_nonzero_poly(numerator(s))
    else
      print("\\frac{")
      _factor_and_print_nonzero_poly(numerator(s))
      print("}{")
      _factor_and_print_nonzero_poly(denominator(s))
      print("}")
    end
    return
  end
end

function _factor_and_print_nonzero_poly(p)
  factored = factor(p)
  u = unit(factored)
  if u != 1
    if u == -1
      print("-")
    else
      print("($u)")
    end
  end
  for f in factored
    print("(")
    print(_format_polynomial("$(f[1])"))
    print(")")
    if f[2] > 1
      print("^{$(f[2])}")
    end
  end
end

function _format_polynomial(poly::String)
  return replace(poly, r"\bt(\d+)\b" => s"t_{\1}")
end

@doc raw"""
    quantum_product_at_q1(G::AbstractGKM_graph, class)

Return the matrix of equivariant quantum multiplication on `G` by the class `class` after setting $q=1$.
All structure constants are computed with respect to the [standard basis](#The-standard-basis) $(e_i)_{i=1}^N$.
Thus, the output is the $N\times N$ matrix expressing the linear map

```math
  \rm{class} \ast|_{q=1} \colon H_T^*(X;\mathbb{Q}) \longrightarrow H_T^*(X;\mathbb{Q}),
  \hspace{7mm}
  a\mapsto (\rm{class} \ast a)|_{q=1}.
```


in the basis $(e_i)_{i=1}^N$.

!!! warning
    This requires `is_strictly_nef(G)==true` as otherwise the quantum product might have infinitely many summands, so setting $q=1$ is not well-defined.

# Example

Recall from the example in [`QH_structure_constants`](@ref) that we computed the equivariant quantum products
of $\mathbb{P}^1$ in the fixed point basis $f_1,f_2$.
To relate this to the output of this function, we need to convert to the standard basis $e_1,e_2$ using
```math
e_1 = \frac{f_1}{t_1-t_2},
\hspace{7mm}
e_2 = \frac{f_2}{-t_1+t_2}
```
and set $q=1$.
The following code confirms that the result is as expected.

```jldoctest quantum_product_at_q1
julia> P1 = projective_space(GKM_graph, 1);

julia> quantum_product_at_q1(P1, point_class(P1, 1))
[(t1^2 - 2*t1*t2 + t2^2 + 1)//(t1 - t2)    1//(t1 - t2)]
[                         -1//(t1 - t2)   -1//(t1 - t2)]

julia> (t1, t2) = gens(P1.equivariantCohomology.coeffRing);

julia> quantum_product_at_q1(P1, [t1, t2])
[(t1^2 - t1*t2 + 1)//(t1 - t2)                    1//(t1 - t2)]
[                -1//(t1 - t2)   (t1*t2 - t2^2 - 1)//(t1 - t2)]
```
"""
function quantum_product_at_q1(G::AbstractGKM_graph, class)

  # If class is given as a vector of coefficients, convert it to a free module element over the
  # parent of its elements (which may be the coefficient ring or its fraction field).
  if class isa AbstractVector
    @req !isempty(class) "class must not be empty"
    @req length(class) == n_vertices(G.g) "class must have one entry per vertex of G"
    P = parent(class[1])
    @req all(c -> parent(c) == P, class) "all entries of class must have the same parent"
    class = free_module(P, length(class))(class)
  end

  SC = QH_structure_constants(G; show_progress=false)
  nv = n_vertices(G.g)

  M = zero_matrix(G.equivariantCohomology.coeffRingLocalized, nv, nv)
  g = gens(G.equivariantCohomology.cohomRingLocalized)

  for beta in keys(SC)
    for i in 1:nv
      cg = quantum_product(G, beta, class, g[i])
      M[i,:] += [cg[j] for j in 1:nv]

    end
  end

  return M
end

@doc raw"""
    c1_at_q1(G::AbstractGKM_graph)

The same as [`quantum_product_at_q1`](@ref)`(G, `[`first_chern_class`](@ref)`(G))`.
That is, return the matrix of the equivariant quantum product with $c_1(T_X)$ at $q=1$, expressed in the [standard basis](#The-standard-basis), 
where `G` is the GKM graph of the GKM space $X$.

# Example

Let us compute the matrix of the equivariant quantum product with $c_1(T_{\mathbb{P}^1})$ on $\mathbb{P}^1$ at $q=1$.

```jldoctest c1_at_q1
julia> c1_at_q1(projective_space(GKM_graph, 1))
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]
```
"""
function c1_at_q1(G::AbstractGKM_graph)
  return quantum_product_at_q1(G, first_chern_class(G))
end

@doc raw"""
    conjecture_O_eigenvalues(G::AbstractGKM_graph; printData::Bool=true)

Return the eigenvalues of quantum multiplication by $c_1^T(TX)$, the equivariant first Chern class of the tangent bundle at $q=1, t=0$,
where $t$ are the equivariant parameters.

!!! warning
    This requires `is_strictly_nef(G)==true` as otherwise the quantum product might have infinitely many summands, so setting $q=1$ is not well-defined.
    It also requires `is_compact(G)==true` as otherwise the non-equivariant limit may be undefined.

# Example

Let us see the examples of $\mathbb{P}^1$, $\mathbb{P}^2$, and $\mathbb{P}^3$.

In the first two julia command lines, we demonstrate the following:
Even though the matrix $M$ of multiplication by $c_1(T_{\mathbb{P}^1})$ has rational functions in the equivariant parameters $t$ as entries,
[`conjecture_O_eigenvalues`](@ref) can still compute the eigenvalues at the non-equivariant limit $t=0$ because the characteristic polynomial
of $M$ has coefficients that are polynomials in $t$.
This follows from the existence of *some* $H_T^*(\text{point};\mathbb{Q})$-linear basis for $H_T^*(X;\mathbb{Q})$
and basis-independence of the characteristic polynomial.

```jldoctest conjecture_O_eigenvalues
julia> c1_at_q1(projective_space(GKM_graph, 1))
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]

julia> characteristic_polynomial(ans)
x^2 - t1^2 + 2*t1*t2 - t2^2 - 4

julia> conjecture_O_eigenvalues(projective_space(GKM_graph, 1)) # P^1
Characteristic poly of c1(TX)* at q=1, t=0:
x^2 - 4
2-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a1: -2.00000}

julia> conjecture_O_eigenvalues(projective_space(GKM_graph, 2)) # P^2
Characteristic poly of c1(TX)* at q=1, t=0:
x^3 - 27
3-element Vector{QQBarFieldElem}:
 {a1: 3.00000}
 {a2: -1.50000 + 2.59808*im}
 {a2: -1.50000 - 2.59808*im}

julia> conjecture_O_eigenvalues(projective_space(GKM_graph, 3)) # P^3
Characteristic poly of c1(TX)* at q=1, t=0:
x^4 - 256
4-element Vector{QQBarFieldElem}:
 {a1: 4.00000}
 {a1: -4.00000}
 {a2: 4.00000*im}
 {a2: -4.00000*im}
```
"""
function conjecture_O_eigenvalues(G::AbstractGKM_graph; printData::Bool=true)
  @req is_compact(G) "G needs to be compact as otherwise we cannot specialize to the non-equivariant limit before computing eigenvalues."
  c1Mat = c1_at_q1(G)
  chi = charpoly(c1Mat)
  chi0 = polynomial(QQ, [0])
  z = repeat([0], rank_torus(G))
  for i in 0:(length(chi)-1)
    set_coefficient!(chi0, i, evaluate(coeff(chi, i), z))
  end
  printData && println("Characteristic poly of c1(TX)* at q=1, t=0:\n$chi0")
  return roots(QQBar, chi0)
end

@doc raw"""
    QH_is_associative(G::AbstractGKM_graph; printDiagnostics::Bool) -> Bool

Return whether the calculated structure constants of $QH_T^*(X)$ are associative.
If `G` is the GKM graph of a GKM variety or Hamiltonian GKM space (see [Definition](../GKM/GKM.md#Definition)), then this should always return `true`.

!!! warning
    This requires `is_strictly_nef(G)==true` as otherwise there might be infinitely many structure constants to check.

# Optional arguments
 - `printDiagnostics::Bool`: If this is `true` and the function's output is `false`, then the indices where associativity fails are printed.

# Example
```jldoctest QH_is_associative
julia> P3 = projective_space(GKM_graph, 3);

julia> QH_is_associative(P3)
true
```
"""
function QH_is_associative(G::AbstractGKM_graph; printDiagnostics::Bool = true)::Bool #TODO: make this usable with non-nef spaces!
  nv = n_vertices(G.g)
  g = [QH_class(G, x) for x in gens(G.equivariantCohomology.cohomRingLocalized)]
  for i in 1:nv, j in 1:nv, k in 1:nv
    if (g[i]*g[j])*g[k] != g[i]*(g[j]*g[k])
      printDiagnostics && println("Not associative for: ($i, $j, $k)")
      return false
    end
  end
  return true
end

@doc raw"""
    QH_is_commutative(G::AbstractGKM_graph) -> Bool

Return whether the calculated structure constants of $QH_T^*(X)$ are commutative.
If `G` is the GKM graph of a GKM variety or Hamiltonian GKM space (see [Definition](../GKM/GKM.md#Definition)), then this should always return `true`.

!!! warning
    This requires `is_strictly_nef(G)==true` as otherwise there might be infinitely many structure constants to check.

# Example
```jldoctest QH_is_commutative
julia> P3 = projective_space(GKM_graph, 3);

julia> QH_is_commutative(P3)
true
```
"""
function QH_is_commutative(G::AbstractGKM_graph)::Bool #TODO: make this usable with non-nef spaces!
  nv = n_vertices(G.g)
  g = [QH_class(G, x) for x in gens(G.equivariantCohomology.cohomRingLocalized)]
  for i in 1:nv, j in 1:nv
    if g[i] * g[j] != g[j] * g[i]
      return false
    end
  end
  return true
end


@doc raw"""
    QH_is_polynomial(G::AbstractGKM_graph) -> Bool

Return whether all structure constants of the equivariant quantum product of `G` calculated so far are polynomial (rather than fractions of polynomials).
!!! note
    This does not calculate any structure constants afresh but checks all constants calculated so far.
    To calculate them, use [`QH_structure_constants`](@ref) (see above).

# Example
```jldoctest QH_is_polynomial
julia> P3 = projective_space(GKM_graph, 3);

julia> QH_structure_constants(P3; show_progress=false);

julia> QH_is_polynomial(P3)
true
```
"""
function QH_is_polynomial(G::AbstractGKM_graph)::Bool
  S = G.QH_structure_consts
  for b in keys(S)
    for k in keys(S[b])
      s = S[b][k]
      if !_is_polynomial(s)
        println("QH is not polynomial for curve class $b at index $k.")
        return false
      end
    end
  end
  return true
end

"""
    QH_is_homogeneous(G::AbstractGKM_graph) -> Bool

Return whether all structure constants of the equivariant quantum product of `G` calculated so far are homogeneous.
!!! note
    This does not calculate any structure constants afresh but checks all constants calculated so far.
    To calculate them, use [`QH_structure_constants`](@ref) (see above).

# Example
```jldoctest
julia> P3 = projective_space(GKM_graph, 3);

julia> QH_structure_constants(P3; show_progress=false);

julia> QH_is_homogeneous(P3)
true
```
"""
function QH_is_homogeneous(G::AbstractGKM_graph)::Bool
  S = G.QH_structure_consts
  for b in keys(S)
    for k in keys(S[b])
      s = S[b][k]
      if !_is_homogeneous(s)
        println("QH is not homogeneous for curve class $b at index $k.")
        return false
      end
    end
  end
  return true
end

"""
    QH_supporting_curve_classes(G::AbstractGKM_graph)

Return a list of all curve classes of `G` in which a non-zero structure constant for the equivariant quantum product has been calculated.
!!! note
    This does not calculate any structure constants afresh but works with all constants calculated so far.
    To calculate them, use [`QH_structure_constants`](@ref).
"""
function QH_supporting_curve_classes(G::AbstractGKM_graph)
  S = G.QH_structure_consts
  res = CurveClass_type[]
  for b in keys(S)
    if !all(s -> iszero(s), S[b])
      push!(res, b)
    end
  end
  return res
end