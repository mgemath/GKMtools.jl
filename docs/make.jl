using Documenter, DocumenterCitations
using GKMtools


bib = CitationBibliography(
    joinpath(@__DIR__, "src", "gkm_references.bib");
    style=:alpha
)

DocMeta.setdocmeta!(GKMtools, :DocTestSetup, :(using Oscar, GKMtools); recursive=true)


pages = [
        "Home" => "index.md",
        "GKM spaces" => ["GKM Graphs" => "GKM/GKM.md", 
                            "Constructors" => "GKM/Constructors.md", 
                            "Properties" => "GKM/Properties.md",
                            "Connections" => "GKM/Connections.md", 
                            "Standard Constructions" => "GKM/STDconstructions.md",
                            "Low dimensional examples" => "GKM/LowdimExamples.md",
                            "Operators" => "GKM/Operators.md", 
                            "Cohomology" => "GKM/Cohomology.md",
                            "Curve Classes" => "GKM/CurveClasses.md",
                            "Vector Bundles" => "GKM/Vectorbundles.md",
                            "Seidel Space" => "GKM/Seidelspace.md"],
        "Gromov--Witten theory & Quantum Cohomology" => ["Gromov--Witten invariants" => "GW/GW.md",
                                                "Quantum Cohomology" => "GW/QH.md",
                                                "Seidel Elements / Shift Operators" => "GW/SeidelElements.md"],
        "Examples from the Article" => [
            "BPS Numbers (Table 5.1)" => "Article/BPS.md",
            "Calabi-Yau Edges (Figure 5.2)" => "Article/Fig5_2.md",
            "Twisted Flag Manifold (Figure 5.3 and Theorem 5.13)" => "Article/Twisted_flag.md",
            "Realizability (Section 5.3.1)" => "Article/Realizability.md"
        ],
        "Miscellaneous" => "Misc/Misc.md",
        "References" => "references.md"]

makedocs(
    sitename = "GKMtools",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1
    ),
    modules = [GKMtools],
    warnonly = true,
    pages = pages,
    plugins = [bib],
    doctest = false,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mgemath/GKMtools.jl.git",
    devbranch = "master",
    devurl = "dev",      # the URL for the development branch
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)
