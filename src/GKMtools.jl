@doc raw"""
**GKMtools** is a Julia package for GKM spaces and their equivariant Gromov–Witten theory.

"""

module GKMtools

using Oscar, Combinatorics, ProgressMeter, nauty_jll, Graphs
using Serialization # optional, for the possibility to store data in the future

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

# Load imports from Oscar
include("imports.jl")

# Export functions
include("exports.jl")

# Types
include("types.jl")

# Core combinatorial data and graph construction
include("core/core.jl")

# Homology
include("homology/homology.jl")

# Cohomology
include("cohomology/cohomology.jl")

# Quantum
include("quantum/quantum.jl")

# Connection
include("connection/connection.jl")

# Smooth
include("smooth/smooth.jl")

# Orbifold
include("orbifold/orbifold.jl")

# Operators
include("operators/operators.jl")

# GKM subgraphs and blowups
include("subgraph/subgraph.jl")

# Vector bundles
include("vb/vb.jl")

# Examples
include("examples/examples.jl")

# Gromov Witten
include("GW/GW.jl")

end # module
