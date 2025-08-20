# GKMtools

A Julia package for computations in GKM theory
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/GKMtools.jl/stable/) [![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://mgemath.github.io/GKMtools.jl/dev/)

## Installation
This package requires Oscar, so make sure that you can use Oscar before installing this package. See https://www.oscar-system.org/install/.
In order to install the most stable version of this package, type:
```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="v0.11.1")
```
For the version under development, type:
```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/mgemath/GKMtools.jl", rev="master")
```
After the installation, simply type:
```julia-repl
julia> using Oscar, GKMtools
```
every time you want to use the program.