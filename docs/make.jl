using NetworkHistogram
using Documenter

DocMeta.setdocmeta!(NetworkHistogram, :DocTestSetup, :(using NetworkHistogram); recursive=true)

makedocs(;
    modules=[NetworkHistogram],
    authors="Jake Grainger, Charles Dufour",
    repo="https://github.com/SDS-EPFL/NetworkHistogram.jl/blob/{commit}{path}#{line}",
    sitename="NetworkHistogram.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SDS-EPFL.github.io/NetworkHistogram.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SDS-EPFL/NetworkHistogram.jl",
    devbranch="main",
)
