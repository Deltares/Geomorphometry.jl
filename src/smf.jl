"""
```
B = smf(A; ω, slope, dhₘ, dh₀, cellsize)
```
Applies the simple morphological filter by [Pingel et al. (2013)](@cite pingelImprovedSimpleMorphological2013a) to `A`.

# Output
- `B::Array{Float64,2}` A filtered version of A

# Arguments
- `A::Array{T,2}` Input Array
- `ω::Float64=18.` Maximum window size [m]
- `slope::Float64=0.01` Terrain slope [m/m]
- `cellsize::Float64=1.` Cellsize in [m]
"""
function smf(
    A::AbstractMatrix{<:Real};
    ω::Real = 17.0,
    slope::Real = 0.01,
    cellsize = abs(first(cellsize(A))),
)
    out = similar(A, Union{Missing, nonmissingtype(eltype(A))})
    out .= A
    lastsurface = -copy(A)
    is_low = falses(size(A))
    radii = 1:(round(Int, ω / cellsize) >> 1)

    threshold = slope * cellsize
    thissurface = opening_circ(lastsurface, 1)
    is_low .|= ((lastsurface - thissurface) .> threshold)
    lastsurface .= thissurface

    lastsurface = copy(A)
    is_obj = falses(size(A))
    for radius in radii
        threshold = slope * radius * cellsize
        thissurface = opening_circ(lastsurface, radius)
        is_obj .|= ((lastsurface - thissurface) .> threshold)
        lastsurface .= thissurface
    end
    out[is_low .| is_obj] .= missing
    return out
end
