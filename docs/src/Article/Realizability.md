# Realizability (Example 5.3.1)

We provide here the code producing the calculations in Section 5.3.1 of the accompanying article.

```jldoctest realizability_Ex_5_3_1
julia> G = gkm_2d([1 0; 0 1; -1 0; 0 -1; 1 0; 0 1; -1 0; 0 -1;])
GKM graph with 8 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (0, -1)
4 -> 3 => (1, 0)
5 -> 4 => (0, 1)
6 -> 5 => (-1, 0)
7 -> 6 => (0, -1)
8 -> 1 => (0, -1)
8 -> 7 => (1, 0)

julia> beta = curve_class(G, Edge(1, 2)) + curve_class(G, Edge(2, 3))
(1, 1, -1, -1, 1, 1)

julia> gromov_witten(G, beta, 0, class_one())
1//(t1^2*t2 - t1*t2^2)

julia> p1 = point_class(G, 1)
t1*t2*e[1]

julia> quantum_product(G, beta, p1, p1)
(t2^2//(t1^2 - t1*t2), -t2//t1, -t2//(t1 - t2), 0, 0, 0, 0, 0)
```