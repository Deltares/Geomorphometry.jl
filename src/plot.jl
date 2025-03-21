# using RecipesBase
"""
    image = pssm(dem; exaggeration=2.3, resolution=1.0)

Perceptually Shaded Slope Map by [Pingel, Clarke., (2014)](@cite pingelPerceptuallyShadedSlope2014a).

# Output
- `image::Gray{T,2}` Grayscale image

# Arguments
- `A::Array{Real,2}` Input Array
- `exaggeration::Real=2.3` Factor to exaggerate elevation
- `cellsize::Real=1.0` Size of cell to account for horizontal resolution if different from vertical resolution
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
    hillshade(dem::Matrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=cellsize(dem))

hillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref) given a light source with azimuth and zenith angles in °, as defined in
[Burrough, P. A., and McDonell, R. A., (1998)](@cite burroughPrinciplesGeographicalInformation2015).
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
    multihillshade(dem::AbstractMatrix{<:Real}; cellsize=cellsize(dem))

multihillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref). Like [`hillshade`](@ref), but now using multiple sources as defined in
[Mark, R.K. (1992)](@cite mark1992multidirectional), similar to GDALs -multidirectional.
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
