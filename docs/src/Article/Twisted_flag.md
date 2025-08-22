# Twisted Flag Manifold (Figure 5.3 and Theorem 5.13)

In this section, we provide the code that generates $H_2(X;\mathbb{Z})$ and $\int_{C_e} c_1(T_X)$ for Figure 5.3 in the article, where $X$ is the twisted flag manifold.
We also probide

## Figure 5.3

Conveniently, the GKM graph of the twisted flag manifold is implemented as a single commend.
Hence, we obtain Figure 5.3 using the following code:

```jldoctest Twisted_flag_figure_5_3
julia> G = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)

julia> print_curve_classes(G)
2 -> 1: (0, 1), Chern number: 4
3 -> 2: (-1, 1), Chern number: 2
4 -> 1: (1, 0), Chern number: 2
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 2
5 -> 4: (-1, 1), Chern number: 2
6 -> 1: (1, 1), Chern number: 6
6 -> 3: (1, 0), Chern number: 2
6 -> 5: (0, 1), Chern number: 4
```

The curve classes $\beta$ and $\gamma$ in Figure 5.2 correspond to $(-2, 1)$ and $(1, 0)$, respectively, in the output.

## Theorem 5.13

Continuing from the setup above, the following code produces the integer coefficients in (5.14) in Theorem 5.13.

```julia
b = curve_class(G, Edge(3, 4))
g = curve_class(G, Edge(1, 4))

p = point_class(G, 1)

for d in 0:5
  gw = quantum_product(G, 3*g+d*b, p, p; useStructureConstants=false, fastMode=true, distantVertex = 5)
  println("Coefficient of q^{2\\beta + $d*\\gamma}: $gw")
end
```

Note that the final three optional arguments are merely to optimize performance for this particular GKM graph.