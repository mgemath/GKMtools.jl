# Enumeration of curves contained in a hyperplane (Section 5.4)

We provide here the code of the calculations described in Section 5.4 of the accompanying article. The problem we want to solve is the enumeration of curves in the projective space $\mathbb{P}^n$ meeting linear subspaces and contained in a $k$-dimensional linear subspace.

The incident variety $S\subsetneq \mathbb{P}^n \times G(n, n+1)$ is the flag variety $F(1,k,n-k)$. See [Standard Constructions](../GKM/STDconstructions.md) for our convention on the notation of flag varieties. Let us consider the case $n=3$ and $k=2$.

```jldoctest Planar_c; setup = :(using Oscar, GKMtools) 
julia> n = 3;

julia> k = 2;

julia> S = flag_variety(GKM_graph, [1, k, n-k]);
```
Now, let us define the pull-backs of a point and of a line in $\mathbb{P}^3$.

```jldoctest Planar_c
julia> vert = D -> ["$a$b$c" for (a,b,c) in Iterators.product(Iterators.repeated(1:4, 3)...) if (a<=D && b < c && b != a && c != a)];

julia> pb_point = poincare_dual(gkm_subgraph_from_vertices(S, vert(1)));

julia> pb_line = poincare_dual(gkm_subgraph_from_vertices(S, vert(2)));
```
We want to compute the number of curves of degree $d$ in $\mathbb{P}^3$ meeting $r$ lines and $s$ points, assumed in general position. We must have $m=r+s$ and $r+2s=3d+2$, where $m$ is the number of marks. Let us compute the case of cubics, that is $d=3$.

```jldoctest Planar_c
julia> s, r = 5, 1; # planar cubics meeting 5 points and 1 line

julia> m = r + s;

julia> beta = curve_class(S, "123", "213"); # curve class of a line

julia> d = 3; # we are taking cubics

julia> P = prod(i -> ev(i, pb_line), 1:r; init = 1) * prod(i -> ev(i, pb_point), (r+1):m; init = 1);

julia> gromov_witten(S, d * beta, m, P; fast_mode = true, show_bar = false)
0
```

Note that this result is expected, as there are no planar cubics meeting $5$ points. In order to get non-trivial results, we can change the values of $r$ and $s$. For example, 

```julia-repl
julia> s, r = 1, 9; # planar cubics meeting 1 point and 9 lines

julia> m = r + s;

julia> P = prod(i -> ev(i, pb_line), 1:r; init = 1) * prod(i -> ev(i, pb_point), (r+1):m; init = 1);

julia> gromov_witten(S, d * beta, m, P; fast_mode = true) # about an hour in our machine
1392
```

We may consider the same problem with $n=4$ and $k=2$. That is, planar cubics in $\mathbb{P}^4$ meeting $s$ points, $r$ lines and $p$ planes. Here an example:

```julia-repl
julia> n = 4;

julia> k = 2;

julia> S = flag_variety(GKM_graph, [1, k, n-k]);

julia> beta = curve_class(S, "123", "213");

julia> vert = D -> ["$a$b$c" for (a,b,c) in Iterators.product(Iterators.repeated(1:5, 3)...) if (a<=D && b < c && b != a && c != a)];

julia> pb_point = poincare_dual(gkm_subgraph_from_vertices(S, vert(1)));

julia> pb_line = poincare_dual(gkm_subgraph_from_vertices(S, vert(2)));

julia> pb_plane = poincare_dual(gkm_subgraph_from_vertices(S, vert(3)));

julia> s, r, p = 0, 6, 2; #planar cubics meeting 6 lines and 2 planes

julia> m = s + r + p;

julia> P = prod(i -> ev(i, pb_point), 1:s; init=1) * prod(i -> ev(i, pb_line), (s+1):(s+r);init=1) * prod(i -> ev(i, pb_plane), (s+r+1):m;init=1);

julia> d = 3;

julia> gromov_witten(S, d * beta, m, P; fast_mode = true, show_bar = true) # about 17 minutes in our machine
60
```

On the other hand, if we change the last lines, we get the following,

```julia-repl
julia> s, r, p = 0, 5, 4; #planar cubics meeting 5 lines and 4 planes

julia> m = s + r + p;

julia> P = prod(i -> ev(i, pb_point), 1:s; init=1) * prod(i -> ev(i, pb_line), (s+1):(s+r);init=1) * prod(i -> ev(i, pb_plane), (s+r+1):m;init=1);

julia> gromov_witten(S, d * beta, m, P; fast_mode = true, show_bar = true)  # about an hour in our machine
540
```