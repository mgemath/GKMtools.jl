# GKMtools.jl
*A Julia package for GKM spaces and their equivariant Gromov–Witten theory.*

This Julia package is work in progress. It offers support for calculations involving GKM spaces, including their equivariant Gromov–Witten invariants.

It comes with the supporting article: Daniel Holmes and Giosuè Muratore, Computations in Equivariant Gromov–Witten theory of GKM spaces.

The package is divided in two parts. The first one deals with foundational material on [GKM spaces](GKM/GKM.md) in general. The second part is dedicated to the computation of [equivariant Gromov–Witten invariants](GW/GW.md), [equivariant quantum cohomology](GW/QH.md) and [equivariant Seidel elements (shift operators)](GW/SeidelElements.md).

We also include in this documentation all [examples from the article](Article/BPS.md) where this is necessary to make them reproducible.

## Installation
This package depends on **Oscar**, so make sure that Oscar is installed and runs correctly on your system (see the [installation instructions](https://www.oscar-system.org/install/)). **Important:** GKMtools.jl must be installed in the same environment where Oscar is available. For example, if Oscar is installed inside WSL on Microsoft Windows, then GKMtools.jl must also be installed and run within that same WSL distribution.

To install the latest stable release of this package, run:

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="v0.13.0")
```

To install the development version (a preview of the upcoming **v1.0.0**), use:

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="master")
```

Once installed, load the package alongside Oscar with:

```julia-repl
julia> using Oscar, GKMtools
```

Copyright (c) 2025: [Daniel Holmes](https://www.daniel-holmes.at/) and [Giosuè Muratore](https://sites.google.com/view/giosue-muratore)