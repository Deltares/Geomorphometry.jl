module GeoArrayOps

using ImageFiltering

include("utils.jl")
include("pmf.jl")
include("smf.jl")
include("plot.jl")
include("spread.jl")

export pmf
export smf
export pssm
export spread
export spread2

precompile(pmf, (Matrix{Float64},))
precompile(smf, (Matrix{Float64},))
precompile(pssm, (Matrix{Float64},))
precompile(spread, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))
precompile(spread2, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))
precompile(spread, (Matrix{Float64}, Matrix{Float64}, Float64))
precompile(spread, (Matrix{Float64}, Float64, Float64))

end # module
