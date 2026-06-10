"""
```
B = simple_morphological_filter(A; ω, slope, cellsize)
```
Applies the simple morphological filter by [Pingel et al. (2013)](@cite pingelImprovedSimpleMorphological2013a) to `A`.

# Output
- `B::Matrix{Union{Missing,T}}` A copy of `A` with filtered (non-ground) cells set to `missing`.

# Arguments
- `A::AbstractMatrix{<:Real}` Input Array
- `ω::Real=17.` Maximum window size [m]
- `slope::Real=0.01` Terrain slope [m/m]
- `cellsize=abs(first(cellsize(A)))` Cellsize in [m]
"""
function simple_morphological_filter(
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
@deprecate smf simple_morphological_filter
