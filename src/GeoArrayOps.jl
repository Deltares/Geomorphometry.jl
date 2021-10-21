module GeoArrayOps

using ImageFiltering

include("utils.jl")
include("pmf.jl")
include("smf.jl")
include("plot.jl")

export pmf
export smf
export pssm

precompile(pmf, (Matrix{Float64},))
precompile(smf, (Matrix{Float64},))
precompile(pssm, (Matrix{Float64},))

end # module
