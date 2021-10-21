module GeoRasterFiltering

using ImageFiltering

include("utils.jl")
include("pmf.jl")
include("smf.jl")

export pmf
export smf

end # module
