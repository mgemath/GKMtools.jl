# GKMtools

A Julia package for computations in GKM theory
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/GKMtools.jl/stable/) [![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://mgemath.github.io/GKMtools.jl/dev/)

This Julia package is work in progress. It offers support for calculations involving GKM spaces, including their equivariant Gromov–Witten invariants.

It comes with the supporting article: Daniel Holmes and Giosuè Muratore, Computations in Equivariant Gromov–Witten theory of GKM spaces.

The package is divided in two parts. The first one deals with foundational material on [GKM spaces](GKM/GKM.md) in general. The second part is dedicated to the computation of [equivariant Gromov–Witten invariants](GW/GW.md), [equivariant quantum cohomology](GW/QH.md) and [equivariant Seidel elements (shift operators)](GW/SeidelElements.md).

We also include in this documentation all [examples from the article](Article/BPS.md) where this is necessary to make them reproducible.

We refer to the documentation page for [installation instructions](https://mgemath.github.io/GKMtools.jl/stable/). 

## Current state of the `flags` branch

On the branch `flags`, we gradually introduce the feature of GKM graphs with standalone flags that do not belong to some edge.
We record here the status of this change.

#### Files with full flag support:

Checked by Daniel:

- `Types.jl`
- `GKMgraphs.jl`
- `GKMconnections.jl`
- `GKMsubgraphs.jl`
- `cohomology.jl`

Not yet checked (i.e. next TODOs for Daniel):

- `product.jl` D: check!
- `equivariant_bundles.jl` D: check, in particular `total_space` and `projectivization`.

#### Files with partial flag support:

- `blowup.jl`: 
    - construction should work (D: check!)
    - connection needs to be re-implemented.

#### Not yet updated:

Definitely broken and needs update:

- `Seidel_space.jl` uses old connection generation

Probably still works for the compact case, but needs update:

- `betti.jl`
- `different_w_types.jl`
- Everything in `GW/`

#### Leftovers for backwards compatibility:

These may be removed later:
- `Types.jl`
    - The field `w` is redundant with `weights_at_vertex` and should be removed eventually.

#### Docs to be updated:
The entire docs need a readthrough and update to reflect the mathematical level of generality we have.
Some particular points:

- `Connections.md` Introduction needs update for flags

#### Files that work unchanged:

All other Julia files should work unchanged, but thorough testing is required. Doctests will help after human checks.