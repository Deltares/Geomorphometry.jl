__precompile__()
module GeoArrayOps
using ProgressMeter

include("utils.jl")
include("pmf.jl")
include("smf.jl")
include("psf.jl")
include("plot.jl")
include("spread.jl")
include("terrain.jl")

export pmf, smf, psf
export pssm
export pitremoval
export spread, spread2
export roughness, TRI, TPI

# function __init__()
#     A = rand(Float64, 3, 3)
#     pssm(A)
#     spread(A, A, A)
#     spread2(A, A, A)
#     spread(A, A, 1.)
#     spread(A, 1., 1.)
#     roughness(A)
#     TRI(A)
#     TPI(A)
# end

precompile(pmf, (Matrix{Float64},))
precompile(psf, (Matrix{Float64},))
precompile(pitremoval, (Matrix{Float64},))
precompile(smf, (Matrix{Float64},))
precompile(pssm, (Matrix{Float64},))
precompile(spread, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))
precompile(spread2, (Matrix{Float64}, Matrix{Float64}, Matrix{Float64}))
precompile(spread, (Matrix{Float64}, Matrix{Float64}, Float64))
precompile(spread, (Matrix{Float64}, Float64, Float64))
precompile(roughness, (Matrix{Float64},))
precompile(TRI, (Matrix{Float64},))
precompile(TPI, (Matrix{Float64},))

end # module
