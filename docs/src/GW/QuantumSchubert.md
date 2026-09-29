# Quantum Schubert calculus

For homogeneous spaces, quantum Chevalley allows to compute **ordinary small quantum**
Schubert structure constants from root data, without stable-map localization.
Here we implemented the formula in order to reconstruct all ordinary quantum products.

```julia-repl
julia> R = root_system(:A, 6);

julia> S = [2, 3, 4, 5, 6]; # Correspond to the variety P^6

julia> Q = quantum_schubert_context(R, S); # this stores only the data needed to the Quantum Chevalley

julia> length(Q.representatives)  # number of vertices
7

julia> Q.novikov_degrees # one component of at most degree 7
1-element Vector{Int64}:
 7

julia> Q.representatives # minimal length elements of the Weyl group of G/P
7-element Vector{WeylGroupElem}:
 id
 s1
 s2 * s1
 s3 * s2 * s1
 s4 * s3 * s2 * s1
 s5 * s4 * s3 * s2 * s1
 s6 * s5 * s4 * s3 * s2 * s1

julia> Q.codimensions # the codimention of each representatives
7-element Vector{Int64}:
 0
 1
 2
 3
 4
 5
 6

julia> Q.representatives[1]
id

julia> Q.representatives[2]
s1
```

Let us compute the quantum product `id * s1`. The expected result is `s1`
```julia-repl
julia> quantum_schubert_product(Q, 1, 2)
Dict{Tuple{Int64, Tuple{Vararg{Int64}}}, QQFieldElem} with 1 entry:
  (2, (0,)) => 1
```
The result is a dictionary with entries of the form: `(vertex, degree) => c`. It means the following: The quantum product is the sum of the entries in the form: $c\cdot \sigma_{vertex}\cdot q^{degree}$, where $\sigma_{vertex}$ is the Schubert class given by the vertex. Thus $s_1 = 1 \sigma_{2} q^0$ as expected. Let us try to compute $s_1^9=\sigma_3 * \sigma_7=q=1\cdot\sigma_2\cdot q^1$:
```julia-repl
julia> quantum_schubert_product(Q, 3, 7)
Dict{Tuple{Int64, Tuple{Vararg{Int64}}}, QQFieldElem} with 1 entry:
  (2, (1,)) => 1

julia> quantum_schubert_coefficient(Q, 1, 7, 2, 1)
1
```
We can compute the entire multiplication table of a Schubert class, and store them in a file for later use.
```julia-repl
julia> v = 2;

julia> products = all_quantum_schubert_products(Q, v) # compute all products of v
7-element Vector{Dict{Tuple{Int64, Tuple{Vararg{Int64}}}, BigInt}}:
 Dict((2, (0,)) => 1)
 Dict((3, (0,)) => 1)
 Dict((4, (0,)) => 1)
 Dict((5, (0,)) => 1)
 Dict((6, (0,)) => 1)
 Dict((7, (0,)) => 1)
 Dict((1, (1,)) => 1)

julia> products = serialize_quantum_schubert_products("filename.jls", Q, 2); # compute the product and store in a file
```
We can have it in form of a matrix, use the following. Use `at_q1=true` to impose `q=1`.

```julia-repl
julia> M = quantum_schubert_matrix(products;  at_q1=false)
[0   0   0   0   0   0   q1]
[1   0   0   0   0   0    0]
[0   1   0   0   0   0    0]
[0   0   1   0   0   0    0]
[0   0   0   1   0   0    0]
[0   0   0   0   1   0    0]
[0   0   0   0   0   1    0]

julia> M = quantum_schubert_matrix(products;  at_q1=true)
[0   0   0   0   0   0   1]
[1   0   0   0   0   0   0]
[0   1   0   0   0   0   0]
[0   0   1   0   0   0   0]
[0   0   0   1   0   0   0]
[0   0   0   0   1   0   0]
[0   0   0   0   0   1   0]
```

## API

```@docs
QuantumSchubertContext
quantum_schubert_context
quantum_chevalley_product
quantum_schubert_coefficient
quantum_schubert_product
quantum_schubert_table
serialize_quantum_schubert_products
quantum_schubert_matrix
```

## Products without serialization

```@docs
all_quantum_schubert_products
```
