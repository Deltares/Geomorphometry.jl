# using RecipesBase
"""
    pssm(dem::AbstractMatrix{<:Real}; exaggeration=2.3, cellsize=cellsize(dem), method=Horn())

Perceptually Shaded Slope Map by [Pingel and Clarke (2014)](@cite pingelPerceptuallyShadedSlope2014a).

# Output
- `Matrix{Float32}` Exaggerated [`slope`](@ref) values, suitable for visualization as a grayscale image.

# Arguments
- `dem::AbstractMatrix{<:Real}` Input digital elevation model
- `exaggeration::Real=2.3` Factor to exaggerate elevation
- `cellsize=cellsize(dem)` Cell size, to account for horizontal resolution if different from vertical resolution
- `method::DerivativeMethod=Horn()` Derivative estimator used for the underlying [`slope`](@ref)
"""
function pssm(
    dem::AbstractMatrix{<:Real};
    exaggeration = 2.3,
    cellsize = cellsize(dem),
    method = Horn(),
)
    slope(dem; cellsize, method, exaggeration)
end

"""
    hillshade(dem::AbstractMatrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=cellsize(dem))

hillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref) given a light source with `azimuth` and `zenith` angles in degrees, as defined in
[Burrough, P. A., and McDonell, R. A., (1998)](@cite burroughPrinciplesGeographicalInformation2015).
Returns a `Matrix{Union{Missing,UInt8}}` of illumination values in `0:255`.
"""
function hillshade(
    dem::AbstractMatrix{<:Real};
    azimuth = 315.0,
    zenith = 45.0,
    cellsize = cellsize(dem),
)
    dst = similar(dem, Union{Missing, UInt8})
    zenithr = deg2rad(zenith)
    azimuthr = deg2rad(azimuth)

    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx, δzδy = (v[3] - v[4]) / (8 * v[5][1]), (v[1] - v[2]) / (8 * v[5][2])
        if δzδx != 0
            a = atan(δzδx, δzδy)
            if a < 0
                a += 2π
            end
        else
            a = π / 2
            if δzδy < 0
                a += 2π
            end
        end
        slope = atan(√(δzδx^2 + δzδy^2))
        something = max(
            0,
            255 * (
                (cos(zenithr) * cos(slope)) +
                (sin(zenithr) * sin(slope) * cos(azimuthr - a))
            ),
        )
        d[i] = isfinite(something) ? round(UInt8, max(0, something)) : missing
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

"""
    multihillshade(dem::AbstractMatrix{<:Real}; azimuth=[225, 270, 315, 360], zenith=45.0, cellsize=cellsize(dem))

multihillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref). Like [`hillshade`](@ref), but combining multiple light sources at the given
`azimuth` angles (degrees) as defined in [Mark, R.K. (1992)](@cite mark1992multidirectional),
similar to GDAL's -multidirectional. Returns a `Matrix{Union{Missing,UInt8}}` of illumination
values in `0:255`.
"""
function multihillshade(
    dem::AbstractMatrix{<:Real};
    azimuth = [225, 270, 315, 360],
    zenith = 45.0,
    cellsize = cellsize(dem),
)
    dst = similar(dem, Union{Missing, UInt8})
    zenithr = deg2rad(zenith)

    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx, δzδy = (v[3] - v[4]) / (8 * v[5][1]), (v[1] - v[2]) / (8 * v[5][2])
        if δzδx != 0
            a = atan(δzδx, δzδy)
            if a < 0
                a += 2π
            end
        else
            a = π / 2
            if δzδy < 0
                a += 2π
            end
        end
        slope = atan(√(δzδx^2 + δzδy^2))

        w225 = sin(a - deg2rad(225 - 90))^2
        w270 = sin(a - deg2rad(270 - 90))^2
        w315 = sin(a - deg2rad(315 - 90))^2
        w360 = sin(a - deg2rad(360 - 90))^2

        α = cos(zenithr) * cos(slope)
        β = sin(zenithr) * sin(slope)
        weights = 0
        something = 0
        for az in azimuth
            weight = sin(a - deg2rad(az - 90))^2
            weights += weight
            something += weight * (α + β * cos(deg2rad(az) - a))
        end
        something /= weights

        d[i] = isfinite(something) ? round(UInt8, max(0, 255 * something)) : missing
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end
