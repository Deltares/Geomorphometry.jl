"""
```
B = smf(A; ω, slope, dhₘ, dh₀, cellsize)
```
Applies the simple morphological filter by *Pingel et al. (2013)* [^pingel2013] to `A`.

# Output
- `B::Array{Float64,2}` A filtered version of A

# Arguments
- `A::Array{T,2}` Input Array
- `ω::Float64=18.` Maximum window size [m]
- `slope::Float64=0.01` Terrain slope [m/m]
- `cellsize::Float64=1.` Cellsize in [m]

[^pingel2013]: Pingel, Thomas J., Keith C. Clarke, and William A. McBride. 2013. ‘An Improved Simple Morphological Filter for the Terrain Classification of Airborne LIDAR Data’. ISPRS Journal of Photogrammetry and Remote Sensing 77 (March): 21–30. <https://doi.org/10.1016/j.isprsjprs.2012.12.002>.
"""
function smf(A::AbstractMatrix{T};
    ω::Real=17.,
    slope::Real=0.01,
    cellsize::Real=1.0) where T <: Real

    lastsurface = -copy(A)
    is_low = falses(size(A))
    radii = 3:2:round(Int, ω / cellsize)
    for radius ∈ radii
        threshold = slope * radius * cellsize
        thissurface = opening_circ(lastsurface, radius)
        is_low .|= ((lastsurface - thissurface) .> threshold)
        lastsurface = thissurface
    end

    lastsurface = copy(A)
    is_obj = falses(size(A))
    for radius ∈ radii
        threshold = slope * radius * cellsize
        thissurface = opening_circ(lastsurface, radius)
        is_obj .|= ((lastsurface - thissurface) .> threshold)
        lastsurface = thissurface
    end
    out = Array{Union{Missing,T}}(A)
    out[is_low .| is_obj] .= missing
    return out
end
