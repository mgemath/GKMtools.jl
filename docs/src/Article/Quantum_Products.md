# Examples of "Effective Gromov-Witten computations for GKM spaces and a semisimple application"

In this section, we reproduce the quantum product computations whose code was printed explicitly in [HM26](@cite).

```jldoctest article2_qprod_ex1; setup = :(using Oscar, GKMtools)
julia> R = root_system(:C, 3);

julia> I = [1, 3];

julia> X = generalized_gkm_flag(R, I);

julia> Q = quantum_schubert_context(X);

julia> products = all_quantum_schubert_products(Q, 2; show_bar=false);

julia> M = quantum_schubert_matrix(products; at_q1 = false)
[0   0   0   0   0   0    0   0   q1   0    0    0]
[1   0   0   0   0   0   q1   0    0   0    0    0]
[0   1   0   0   0   0    0   0    0   0    0    0]
[0   1   0   0   0   0    0   0    0   0   q1    0]
[0   0   1   1   0   0    0   0    0   0    0    0]
[0   0   0   0   2   0    0   1    0   0    0    0]
[0   0   0   0   0   1    0   0    1   0    0    0]
[0   0   0   1   0   0    0   0    0   0    0   q1]
[0   0   0   0   1   0    0   2    0   0    0    0]
[0   0   0   0   0   1    0   0    0   0    0    0]
[0   0   0   0   0   0    1   0    0   1    0    0]
[0   0   0   0   0   0    0   0    0   0    1    0]
```

```jldoctest article2_qprod_ex2; setup = :(using Oscar, GKMtools)
julia> P3 = projective_space(NormalToricVariety, 3);

julia> y = domain(blow_up(P3, 5));

julia> pc = picard_group_with_map(y)[1];

julia> l1 = toric_line_bundle(y, pc([1, 0]));

julia> l2 = toric_line_bundle(y, pc([0, 1]));

julia> L = gkm_vector_bundle_of_toric(l1);

julia> L2 = gkm_vector_bundle_of_toric(l2);

julia> h1 = first_chern_class(L);

julia> h2 = first_chern_class(L2);

julia> Y = baseof(L);

julia> my_basis = [h1^0, h2, h1*h2, h1, h1^2, h1^3];

julia> N = matrix_small_quantum_product(Y, h1; basis = my_basis, show_bar = false)
[0   0   0   0   1   1]
[0   0   0   0   0   1]
[0   1   0   0   0   0]
[1   0   0   0   0   0]
[0   0   0   1   0   0]
[0   0   1   0   1   0]

julia> my_basis = [h1^0, h2, h1*h2, h1, h1^2];

julia> M = twisted_matrix_small_quantum_product(L, h1; basis = my_basis, show_bar = false)
[0      0   1      1   1]
[0      0   1      0   1]
[0   1//2   0   1//2   0]
[1      0   0      0   0]
[0   1//2   0   1//2   0]
```

```julia-repl
julia> V1 = vector_bundle_O(2, [1]);

julia> G, (v1, v2, v3) = (baseof(V1), V1) * (baseof(V1), V1) * (baseof(V1), V1);

julia> V = (v1*v2) + (v1*v3) + (v2*v3);

julia> h1 = chern_class(v1, 1);

julia> h2 = chern_class(v2, 1);

julia> h3 = chern_class(v3, 1);

julia> my_basis = [h1^0, h1, h2, h3, h1^2, h1*h3, h2^2, h1^2*h2];

julia> classes = [h1, 2*h2, 3*h3];

julia> M = twisted_matrix_small_quantum_product(V, classes, basis = my_basis);

julia> discriminant(charpoly(M[1]+M[2]+M[3]))
-3627077673427665895526400000000000000000000000000
```