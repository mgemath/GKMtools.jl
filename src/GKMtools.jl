@doc raw"""
**GKMtest** is a test module for GKM varieties to be included in OSCAR.

"""

module GKMtools

using Oscar, Combinatorics, ProgressMeter

include("imports.jl")
include("exports.jl")

include("GKM/AbstractTypes.jl")
include("GKM/ConcreteTypes.jl")

include("GKM/GKMgraphs.jl")
include("GKM/different_w_types.jl")

include("GKM/betti.jl")

include("GKM/standard_constructions.jl")
include("GKM/GP.jl")

## Low dimensional Examples
include("GKM/lowdimexamples.jl")

include("GKM/Subgraphs.jl")
include("GKM/product.jl")
include("GKM/blowup.jl")


include("GKM/Cohomology.jl")
include("GKM/CurveClasses.jl")
include("GKM/Connections.jl")

## experimental
include("GKM/bruhat.jl")
include("GKM/bott_samelson.jl")

## Miscellaneous
include("misc/bruhatsmoothness.jl")
include("misc/kazhdan_lusztig.jl")

## Seidel Space
include("GKM/Seidel_space.jl")

# ## Vector Bundles
# include("GKM/equivariant_bundles.jl")

###################
## GW invariants ##
###################

## Combinatorial part
include("GW/Combinatorial/Trees.jl")
include("GW/Combinatorial/Colors.jl")
include("GW/Combinatorial/Marks.jl")

## Decorated Trees
include("GW/DecoratedTrees.jl")

## Equivariant Class
include("GW/EquivariantClasses/Rules.jl")
include("GW/EquivariantClasses/Psi_class.jl")
include("GW/EquivariantClasses/class_one_class.jl")
include("GW/EquivariantClasses/ev_class.jl")
include("GW/EquivariantClasses/Euler.jl")

## Gromov-Witten invariants
include("GW/GromovWitten.jl")

end # module GKMtools
