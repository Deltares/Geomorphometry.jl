push!(LOAD_PATH, "../src/")
using Documenter, GeoArrayOps

makedocs(sitename="GeoArrayOps.jl")

deploydocs(
    repo="github.com/Deltares/GeoArrayOps.jl.git",
)
