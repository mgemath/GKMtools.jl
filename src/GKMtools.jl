@doc raw"""
**GKMtools** is a test module for GKM varieties to be included in OSCAR.

"""

module GKMtools

using Oscar, Combinatorics, ProgressMeter, nauty_jll, Graphs

## GKM
include("exports.jl")
include("imports.jl")
include("Types.jl")
include("different_w_types.jl")

## Constructors 
## Properties
include("GKMgraphs.jl")
include("betti.jl")

## Standard Constructions
include("standard_constructions.jl")
include("GP.jl")

## Low dimensional Examples
include("lowdimexamples.jl")

## Operators
include("GKMsubgraphs.jl")
include("product.jl")
include("blowup.jl")

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

## Partial computations
include("GW/partials/partial_comp.jl")

## Positive genus
include("GW/PosGenus/iterator_pos_genus.jl")
include("GW/PosGenus/Main_pos_gen.jl")

end # module GKMtools
