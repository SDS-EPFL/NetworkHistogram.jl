using NetworkHistogram
using Documenter

## Literate preprocessing, maybe move to a separate script later for faster builds

# to run with LiveServer and avoid infinite loops, use
# servedocs(literate_dir=joinpath("docs","literate","tutorials"),skip_dir = joinpath("docs","src","tutorials"))
# adapting `tutorials` to whatever subdir you are working on

using Literate

LITERATE_INPUT = joinpath(@__DIR__, "literate")
LITERATE_OUTPUT = joinpath(@__DIR__, "src")

for dir_path in filter(isdir, readdir(joinpath(@__DIR__, "literate"), join = true))
    dirname = basename(dir_path)

    for (root, _, files) in walkdir(dir_path), file in files
        # ignore non julia files
        splitext(file)[2] == ".jl" || continue
        # full path to a literate script
        ipath = joinpath(root, file)
        # generated output path
        opath = splitdir(replace(ipath, LITERATE_INPUT => LITERATE_OUTPUT))[1]
        # generate the markdown file calling Literate
        Literate.markdown(ipath, opath)
    end
end

DocMeta.setdocmeta!(
    NetworkHistogram, :DocTestSetup, :(using NetworkHistogram);
    recursive = true)

# based on available extensions, include them in the documentation
modules_all = [
    NetworkHistogram,
    Base.get_extension(NetworkHistogram, :BootstrapExt),
    Base.get_extension(NetworkHistogram, :LightMCExt),
    Base.get_extension(NetworkHistogram, :MakieExt),
    Base.get_extension(NetworkHistogram, :PythonOptimalTransport)
]
println("all modules: ", modules_all)

#TODO: safety check, should probably throw an error instead
modules = [filter(!isnothing, modules_all)...]
println("build/docs with modules: ", modules)

using DocumenterInterLinks

links = InterLinks(
    "ot" => "https://pythonot.github.io/",
)

makedocs(;
    modules = modules,
    authors = "Jake Grainger, Charles Dufour",
    repo = "github.com/SDS-EPFL/NetworkHistogram.jl.git",
    sitename = "NetworkHistogram.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://SDS-EPFL.github.io/NetworkHistogram.jl",
        edit_link = "main",
        assets = String[]),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Tutorials" => ["First steps" => "tutorials/simple_graph.md",
            "Multiplex networks" => "tutorials/multiplex_network.md",
            "Weighted networks" => "tutorials/weighted_network.md",
            "Temporal networks" => "tutorials/temporal_networks.md"]],
    checkdocs = :none,
    plugins = [links])

deploydocs(;
    repo = "github.com/SDS-EPFL/NetworkHistogram.jl.git")
