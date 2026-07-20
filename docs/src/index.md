# GKMtools.jl
*A Julia package for GKM spaces and their equivariant Gromov–Witten theory.*

This Julia package is work in progress. It offers support for calculations involving GKM spaces, including their equivariant Gromov–Witten invariants.

It comes with the supporting article: Daniel Holmes and Giosuè Muratore, Computations in Equivariant Gromov–Witten theory of GKM spaces.

The package is divided into two parts. The first one deals with foundational material on [GKM spaces](Generalities/GKM.md), including smooth and orbifold GKM graphs, equivariant cohomology, and equivariant vector bundles. The second part, which is currently being adapted to the new implementation, is dedicated to equivariant Gromov–Witten theory.

## Installation
This package depends on **Oscar**, so make sure that Oscar is installed and runs correctly on your system (see the [installation instructions](https://www.oscar-system.org/install/)). **Important:** GKMtools.jl must be installed in the same environment where Oscar is available. For example, if Oscar is installed inside WSL on Microsoft Windows, then GKMtools.jl must also be installed and run within that same WSL distribution.

To install the latest release of this package, run:

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="master")
```

This package had a substantial remake. To install the latest legacy version (versions **v0.X.X**), use:

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="v0.17.0")
```

Once installed, load the package alongside Oscar with:

```julia-repl
julia> using Oscar, GKMtools
```

Copyright (c) 2025: [Daniel Holmes](https://www.daniel-holmes.at/) and [Giosuè Muratore](https://sites.google.com/view/giosue-muratore)