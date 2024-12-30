# using RecipesBase
"""
    image = pssm(A; exaggeration=2.3, resolution=1.0)

Perceptually Shaded Slope Map by *Pingel, Clarke. 2014* [^pingel2014].

# Output
- `image::Gray{T,2}` Grayscale image

# Arguments
- `A::Array{Real,2}` Input Array
- `exaggeration::Real=2.3` Factor to exaggerate elevation
- `cellsize::Real=1.0` Size of cell to account for horizontal resolution if different from vertical resolution

[^pingel2014]: Pingel, Thomas, and Clarke, Keith. 2014. ‘Perceptually Shaded Slope Maps for the Visualization of Digital Surface Models’. Cartographica: The International Journal for Geographic Information and Geovisualization 49 (4): 225–40. <https://doi.org/10/ggnthv>.
"""
function pssm(
    A::AbstractMatrix{<:Real};
    exaggeration = 2.3,
    cellsize = 1.0,
    method = Horn(),
)
    slope(A; cellsize, method, exaggeration)
end

"""
    hillshade(dem::Matrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=1.0)

hillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref) given a light source with azimuth and zenith angles in °, , as defined in
Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function hillshade(
    dem::AbstractMatrix{<:Real};
    azimuth = 315.0,
    zenith = 45.0,
    cellsize = 1.0,
)
    dst = similar(dem, UInt8)
    zenithr = deg2rad(zenith)
    azimuthr = deg2rad(aspect(azimuth))

    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx, δzδy = (v[1] - v[2]) / (8 * v[5]), (v[3] - v[4]) / (8 * v[5])
        if δzδx != 0
            a = atan(-δzδx, δzδy)
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
        d[i] = round(
            UInt8,
            max(
                0,
                255 * (
                    (cos(zenithr) * cos(slope)) +
                    (sin(zenithr) * sin(slope) * cos(azimuthr - a))
                ),
            ),
        )
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

"""
    multihillshade(dem::Matrix{<:Real}; cellsize=1.0)

multihillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref). Like [`hillshade`](@ref), but now using multiple sources as defined in
https://pubs.usgs.gov/of/1992/of92-422/of92-422.pdf, similar to GDALs -multidirectional.
"""
function multihillshade(dem::AbstractMatrix{<:Real}; cellsize = 1.0)
    dst = similar(dem, UInt8)
    zenithr = deg2rad(60)

    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx, δzδy = (v[1] - v[2]) / (8 * v[5]), (v[3] - v[4]) / (8 * v[5])
        if δzδx != 0
            a = atan(-δzδx, δzδy)
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

        w225 = 0.5 * (1 - cos(2(a - deg2rad(aspect(225)))))
        w270 = 0.5 * (1 - cos(2(a - deg2rad(aspect(270)))))
        w315 = 0.5 * (1 - cos(2(a - deg2rad(aspect(315)))))
        w360 = 0.5 * (1 - cos(2(a - deg2rad(aspect(360)))))

        α = cos(zenithr) * cos(slope)
        β = sin(zenithr) * sin(slope)
        something =
            (
                w225 * (α + β * cos(deg2rad(aspect(225)) - a)) +
                w270 * (α + β * cos(deg2rad(aspect(270)) - a)) +
                w315 * (α + β * cos(deg2rad(aspect(315)) - a)) +
                w360 * (α + β * cos(deg2rad(aspect(360)) - a))
            ) / 2

        d[i] = round(UInt8, max(0, 255 * something))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end
