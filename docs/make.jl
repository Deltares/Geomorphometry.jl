# push!(LOAD_PATH, "../src/")
using Geomorphometry
using Documenter
using DocumenterVitepress
using CairoMakie
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

CairoMakie.activate!(; type = "png")

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
        "Get Started" => "tutorial.md",
        "Concepts" => "concepts.md",
        "Reference" => "reference.md",
        "Guides" => "guide.md",
    ],
    clean = false,
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/Deltares/Geomorphometry.jl.git",
    target = "build",
    devbranch = "master",
    branch = "gh-pages",
    push_preview = true,
)
