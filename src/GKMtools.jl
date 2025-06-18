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

include("GKM/Connections.jl")
# include("GKM/Subgraphs.jl")
# include("GKM/Cohomology.jl")

include("GKM/product.jl")

include("GKM/standard_constructions.jl")
include("GKM/GP.jl")

end # module GKMtools
