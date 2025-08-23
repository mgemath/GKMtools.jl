# Enumeration of curves contained in a hyperplane (Section 5.4)

We provide here the code producing the calculations in Section 5.4 of the accompanying article. The problem we want to solve is the enumeration of curves in the projective space meeting linear subspaces and contained in a hyperplane.

The incident variety $S\subsetneq \mathbb{P}^n \times G(n, n+1)$ is $S=\mathbb{P}(\Omega_{\mathbb{P}^n}(1))\cong \mathbb{P}(\Omega_{\mathbb{P}^n})$. We construct $S$ in the following way:

```jldoctest Planar_c; setup = :(using Oscar, GKMtools) 
julia> P3 = projective_space(GKM_graph, 3);

julia> S = projectivization(cotangent_bd(P3));
GKM graph is valid but not 3-independent, so connections may not be unique.
```
Now, let us define the pull-backs of a point and of a line in $\mathbb{P}^3$.

```jldoctest Planar_c
julia> pb_vertices = dim -> ["[$j]_$i" for i in 1:3 for j in 1:dim];

julia> pb_point = poincare_dual(gkm_subgraph_from_vertices(S, pb_vertices(1)));

julia> pb_line = poincare_dual(gkm_subgraph_from_vertices(S, pb_vertices(2)));
```
If $r$ is the number of lines and $s$ is the number of points, we must have $m=r+s$ and $r+2s=3d+2$, where $m$ is the number of marks and $d$ is the degee of the curves. Let us compute the case of cubics.

```jldoctest Planar_c
julia> s, r = 5, 1; # planar cubics meeting 5 points and 1 line

julia> m = r + s;

julia> beta = curve_class(S, "[2]_2", "[1]_2"); # curve class of a line

julia> d = 3; # we are taking cubics

julia> P = prod(i -> ev(i, pb_line), 1:r) * prod(i -> ev(i, pb_point), (r+1):m);

julia> gromov_witten(S, d * beta, m, P; fast_mode = true, show_bar = false)
0
```

Note that this result is expected, as there are no planar cubics meeting $5$ points. For non trivial results, change the values of $r$ and $s$. For example, 

```julia-repl
julia> s, r = 1, 9; # planar cubics meeting 1 point and 9 lines

julia> m = r + s;

julia> P = prod(i -> ev(i, pb_line), 1:r) * prod(i -> ev(i, pb_point), (r+1):m);

julia> gromov_witten(S, d * beta, m, P; fast_mode = true) # about an hour in our machine
1392
```