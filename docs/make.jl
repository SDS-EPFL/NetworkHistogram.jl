using NetworkHistogram
using Documenter

DocMeta.setdocmeta!(NetworkHistogram, :DocTestSetup, :(using NetworkHistogram);
    recursive = true)

makedocs(;
    modules = [NetworkHistogram],
    authors = "Jake Grainger, Charles Dufour",
    #repo = "github.com/SDS-EPFL/NetworkHistogram.jl.git",
    sitename = "NetworkHistogram.jl",
    #format = Documenter.HTML(;
    #                         prettyurls = get(ENV, "CI", "false") == "true",
    #                         canonical = "https://SDS-EPFL.github.io/NetworkHistogram.jl",
    #                         edit_link = "main",
    #                         assets = String[]),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Optimization hyperparameters" => "rules.md",
        "Development" => "internals.md",
    ],
    checkdocs = :none)

deploydocs(;
    repo = "github.com/SDS-EPFL/NetworkHistogram.jl.git")
