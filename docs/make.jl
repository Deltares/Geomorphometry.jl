push!(LOAD_PATH, "../src/")
using Documenter, GeoArrayOps

makedocs(sitename="GeoArrayOps.jl")

deploydocs(
    repo="github.com/evetion/GeoArrayOps.jl.git",
)
