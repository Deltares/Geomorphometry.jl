module GeoArrayOps

using ImageFiltering

include("utils.jl")
include("pmf.jl")
include("smf.jl")

export pmf
export smf

precompile(pmf, (Matrix{Float64},))
precompile(smf, (Matrix{Float64},))

end # module
