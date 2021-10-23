push!(LOAD_PATH, "../src/")
using Documenter, GeoArrayOps

makedocs(
    sitename="GeoArrayOps.jl",
    pages=[
        "index.md",
        "Tutorials" => [],
        "Topics" => [],
        "Reference" => "reference.md",
        "How-To Guides" => [],
    ]
)

deploydocs(
    repo="github.com/Deltares/GeoArrayOps.jl.git",
)
