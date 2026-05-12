@doc raw"""
**GKMtools** is a Julia package for GKM spaces and their equivariant Gromov–Witten theory.

"""

module GKMtools

using Oscar, Combinatorics, ProgressMeter, nauty_jll, Graphs
using Serialization

const PROJECT_TOML_TEXT = read(joinpath(@__DIR__, "..", "Project.toml"), String)
const VERSION_NUMBER = let
  m = match(r"(?m)^version\s*=\s*\"([^\"]+)\"", PROJECT_TOML_TEXT)
  m === nothing && error("Could not determine the GKMtools version from Project.toml.")
  VersionNumber(only(m.captures))
end
const STARTUP_SUBTITLE = "A Julia package for GKM spaces and their equivariant Gromov-Witten theory."
const DOCS_URL = "https://mgemath.github.io/GKMtools.jl/stable/"
const STARTUP_ART = (
  "   ____ _  ____  __ _              _",
  "  / ___| |/ /  \\/  | |_ ___   ___ | |___",
  " | |  _| ' /| |\\/| | __/ _ \\ / _ \\| / __|",
  " | |_| | . \\| |  | | || (_) | (_) | \\__ \\",
  "  \\____|_|\\_\\_|  |_|\\__\\___/ \\___/|_|___/",
)
const WIDE_DIVIDER = "  --------------------------------------------------"
const AUTHOR_NAMES = "Daniel Holmes and Giosuè Muratore"

_startup_banner_enabled() = lowercase(strip(get(ENV, "GKMTOOLS_STARTUP", "true"))) ∉ ("0", "false", "no", "off")

function _print_banner(io::IO=stdout)
  if displaysize(io)[2] >= 72
    for line in STARTUP_ART
      println(io, line)
    end
    printstyled(io, WIDE_DIVIDER * "\n"; color=:yellow)
    println(io, "  ", STARTUP_SUBTITLE)
    println(io, "  By ", AUTHOR_NAMES)
    println(io, "  Docs: ", DOCS_URL)
    println(io, "  Version: ", VERSION_NUMBER)
    println(io)
  else
    println(io, "GKMtools")
    printstyled(io, WIDE_DIVIDER * "\n"; color=:yellow)
    println(io, STARTUP_SUBTITLE)
    println(io, "By ", AUTHOR_NAMES)
    println(io, "Docs: ", DOCS_URL)
    println(io, "Version: ", VERSION_NUMBER)
    println(io)
  end
end

function __init__()
  isinteractive() || return
  _startup_banner_enabled() || return
  _print_banner()
end

## GKM
include("exports.jl")
include("imports.jl")
include("Types.jl")
include("different_w_types.jl")

## Constructors
## Properties
include("GKMgraphs.jl")
include("betti.jl")
include("indices.jl")

## Standard Constructions
include("standard_constructions.jl")
include("GP.jl")

## Low dimensional Examples
include("lowdimexamples.jl")

## Operators
include("GKMsubgraphs.jl")
include("product.jl")
include("blowup.jl")

## Visualization
include("drawings/admissible_drawings.jl")
include("drawings/convex_drawings.jl")
include("drawings/latex_drawings.jl")

## Connections
include("GKMconnections.jl")

## Cohomology
include("cohomology.jl")
include("curveClasses.jl")

## Vector Bundles
include("equivariant_bundles.jl")

## Seidel Space
include("Seidel_space.jl")


## GW
include("GW/includes.jl")

## obsolate
include("obsolate/obsolate.jl")

## experimental
include("bruhat.jl")
include("bott_samelson.jl")

## Miscellaneous
include("misc/bruhatsmoothness.jl")
include("misc/kazhdan_lusztig.jl")
include("tautological_bd.jl")
include("tautological_bd_GP.jl")
end # module GKMtools
