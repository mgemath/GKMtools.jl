# GKMtools

A Julia package for computations in GKM theory
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/GKMtools.jl/stable/) [![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://mgemath.github.io/GKMtools.jl/dev/)

This Julia package is work in progress. It offers support for calculations involving GKM spaces, including their equivariant Gromov–Witten invariants.

It comes with the supporting article: Daniel Holmes and Giosuè Muratore, Computations in Equivariant Gromov–Witten theory of GKM spaces.

The package is divided in two parts. The first one deals with foundational material on [GKM spaces](GKM/GKM.md) in general. The second part is dedicated to the computation of [equivariant Gromov–Witten invariants](GW/GW.md), [equivariant quantum cohomology](GW/QH.md) and [equivariant Seidel elements (shift operators)](GW/SeidelElements.md).

We also include in this documentation all [examples from the article](Article/BPS.md) where this is necessary to make them reproducible.

We refer to the documentation page for [installation instructions](https://mgemath.github.io/GKMtools.jl/stable/). 

## Current state of the `flags` branch

On the branch `flags`, which is now merged into `GKMtools`, we gradually introduced the feature of GKM graphs with standalone flags that do not belong to some edge.
We record here the status of this change.

#### Files with full flag support:

Checked by Daniel:

- `Types.jl`
- `GKMgraphs.jl`
- `GKMconnections.jl`
- `GKMsubgraphs.jl`
- `cohomology.jl`
- `different_w_types.jl`
- `product.jl`
- `equivariant_bundles.jl`
- `Seidel_space.jl`
- The main function `gromov_witten`.

#### Files with partial flag support:

Checked by Daniel:

- `blowup.jl`: 
    - construction works.
    - natural connection not yet induced, but irrelevant for GW applications.
- `betti.jl`
    - only accepts compact GKM spaces, as combinatorial betti numbers are not well-defined otherwise.
- `equivariant_bundles.jl`
    - `gkm_line_bundle_of_toric` still requires the toric variety to be projective.

#### Not yet updated:

These are features that might just work unchanged, but I am not completely sure.

- Nomarks
- Partials
- Psi classes (also positive genus?)

#### Leftovers for backwards compatibility:

These may be removed later:
- `Types.jl`
    - The field `w` is redundant with `weights_at_vertex` and should be removed eventually.

#### Files that work unchanged:

All other Julia files should work unchanged, but thorough testing is required. All doctests pass at this moment.

#### Other bugfixes or new features:

- Some functions created copies of the GKM graph, but not deep copies.
    - `different_w_types.jl`
    - `substitute_torus` in `GKMgraphs.jl`
The problem was that the resulting fields `curveClasses`, `equivariantCohomology`, and `connection` (and possibly others) carry a reference to their GKM graph, which was still pointing to the old GKM graph.

- `add_edge!` and `add_standalone_flag!` now return the indices (respectively index) of the flag(s) they created. This is useful for certain constructions like `projectivization`.
- `flags_only_gkm_graph` function to create GKM graph with flags but no edges.
- `connect_flags!` to join two standalone flags to an edge in a GKM graph.