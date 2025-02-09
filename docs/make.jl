# push!(LOAD_PATH, "../src/")
using Geomorphometry
using Documenter
using DocumenterVitepress
using CairoMakie
using DocumenterCitations
using Revise
using Downloads

Revise.revise()

# bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))
dir = @__DIR__
# cp(joinpath(dir, "../CHANGELOG.md"), joinpath(dir, "src/CHANGELOG.md"), force = true)
CairoMakie.activate!(; type = "png")

fn = joinpath(dir, "src", "saba.tif")
isfile(fn) && Downloads.download("https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba.tif", fn)

# DocMeta.setdocmeta!(
#     Geomorphometry,
#     :DocTestSetup,
#     :(using Geomorphometry);
#     recursive = true,
# )

makedocs(;
    modules = [Geomorphometry],
    authors = "Maarten Pronk <git@evetion.nl> and contributors",
    repo = "https://github.com/Deltares/Geomorphometry.jl/blob/{commit}{path}#L{line}",
    sitename = "Geomorphometry.jl",
    format = MarkdownVitepress(;
        repo = "github.com/Deltares/Geomorphometry.jl",
        md_output_path = ".",
        build_vitepress = false,
    ),
    doctest = true,
    checkdocs = :all,
    pages = [
        "Home" => "index.md",
        "Getting started" => Any[
            "Installation" => "installation.md",
            "Usage" => "usage.md",
        ],
        "Background" => Any[
            "Concepts" => "concepts.md",
            "Future plans" => "todo.md",
        ],
        "Reference" => Any[
            "API" => "reference.md"
            "Changelog" => "CHANGELOG.md"
        ],
    ],
    clean = false,
    # plugins = [bib],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(;
    repo = "github.com/Deltares/Geomorphometry.jl.git",
    target = "build",
    devbranch = "master",
    branch = "gh-pages",
    push_preview = true,
)
