using Revise
using Geomorphometry
using Documenter
using DocumenterVitepress
using CairoMakie
using DocumenterCitations
using Downloads

Revise.revise()
dir = @__DIR__

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)
cp(joinpath(dir, "../CHANGELOG.md"), joinpath(dir, "src/CHANGELOG.md"); force = true)
CairoMakie.activate!(; type = "png")

fn = joinpath(dir, "src", "saba.tif")
isfile(fn) || Downloads.download(
    "https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba.tif",
    fn,
)
fn = joinpath(dir, "src", "saba_dsm.tif")
isfile(fn) || Downloads.download(
    "https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba_dsm.tif",
    fn,
)
fn = joinpath(dir, "src", "Copernicus_DSM_10_N52_00_E004_00_DEM.tif")
isfile(fn) || Downloads.download(
    "https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/Copernicus_DSM_10_N52_00_E004_00_DEM.tif",
    fn,
)

DocMeta.setdocmeta!(
    Geomorphometry,
    :DocTestSetup,
    :(using Geomorphometry);
    recursive = true,
)

makedocs(;
    modules = [Geomorphometry],
    authors = "Maarten Pronk <git@evetion.nl> and contributors",
    repo = "https://github.com/Deltares/Geomorphometry.jl/blob/{commit}{path}#L{line}",
    sitename = "Geomorphometry.jl",
    format = MarkdownVitepress(;
        repo = "github.com/Deltares/Geomorphometry.jl",
        md_output_path = ".",
        build_vitepress = false,
        devbranch = "main",
    ),
    doctest = true,
    checkdocs = :all,
    pages = [
        "Home" => "index.md",
        "Getting started" => Any[
            "Installation" => "installation.md",
            "Usage" => "usage.md",
            "Experimental" => "experimental.md",
        ],
        "Background" => Any["Concepts" => "concepts.md", "Future plans" => "todo.md"],
        "Reference" => Any[
            "Validation" => "validation.md",
            "API" => "reference.md",
            "Changelog" => "CHANGELOG.md",
            "Bibliography" => "bibliography.md",
        ],
    ],
    clean = false,
    plugins = [bib],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(;
    repo = "github.com/Deltares/Geomorphometry.jl.git",
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
