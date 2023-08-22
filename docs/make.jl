push!(LOAD_PATH, "../src/")
using Documenter, Geomorphometry

makedocs(
    sitename="Geomorphometry.jl",
    pages=[
        "index.md",
        "Tutorials" => [],
        "Topics" => [],
        "Reference" => "reference.md",
        "How-To Guides" => [],
    ]
)

deploydocs(
    repo="github.com/Deltares/Geomorphometry.jl.git",
)
