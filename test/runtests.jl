using Documenter
# include("GKM/GKM.jl")
DocMeta.setdocmeta!(GKMtools, :DocTestSetup, :(using Oscar, GKMtools); recursive=true)
doctest(GKMtools)