using Documenter, GKMtools, Oscar

DocMeta.setdocmeta!(GKMtools, :DocTestSetup, :(using Oscar, GKMtools); recursive=true)
Documenter.doctest(GKMtools; manual=false)
