const neib_8 = @SMatrix [1.0 1 1; 1 0 1; 1 1 1]
const neib_8_inf = @SMatrix [1.0 1 1; 1 Inf 1; 1 1 1]
const neib_8_mask = @SMatrix Bool[1 1 1; 1 0 1; 1 1 1]

const nbkernel = LocalFilters.Kernel{Int8,2}(reshape(1:9, 3, 3))

abstract type DerivativeMethod end

"""Second order finite difference estimator using all 4 neighbors (Zevenbergen and Thorne, 1987)."""
struct ZevenbergenThorne <: DerivativeMethod end

"""Third order finite difference estimator using all 8 neighbors (Horn, 1981)."""
struct Horn <: DerivativeMethod end

"""Maximum Downward Gradient"""
struct MDG <: DerivativeMethod end

"""
    roughness(dem::Matrix{<:Real})

Roughness is the largest inter-cell difference of a central pixel and its surrounding cell, as defined in Wilson et al (2007, Marine Geodesy 30:3-35).
"""
function roughness(dem::AbstractMatrix{<:Real})
    dst = copy(dem)

    initial(a) = (zero(a), a)
    update(v, a, _) = (max(v[1], abs(a - v[2])), v[2])
    store!(d, i, v) = @inbounds d[i] = v[1]

    return localfilter!(dst, dem, 3, initial, update, store!)
end

"""
    TPI(dem::Matrix{<:Real})

TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).
"""
function TPI(dem::AbstractMatrix{<:Real})
    dst = copy(dem)

    initial(a) = zero(a)
    update(v, a, b) = v + a
    store!(d, i, v) = @inbounds d[i] = d[i] - (v - d[i]) / 8

    return localfilter!(dst, dem, 3, initial, update, store!)
end

"""
    TRI(dem::Matrix{<:Real})

TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells.
This algorithm uses the square root of the sum of the square of the difference between a central pixel and its surrounding cells.
This is recommended for terrestrial use cases.
"""
function TRI(dem::AbstractMatrix{<:Real})
    dst = copy(dem)

    initial(a) = (zero(a), a)
    update(v, a, _) = (v[1] + (a - v[2])^2, v[2])
    store!(d, i, v) = @inbounds d[i] = sqrt(v[1])

    return localfilter!(dst, dem, 3, initial, update, store!)
end

"""
    pitremoval(dem::Matrix{<:Real})

Remove pits from a DEM Array if the center cell of a 3x3 patch is `limit` lower or than the surrounding cells.
"""
function pitremoval(dem::AbstractMatrix{<:Real}; limit=0.0)
    dst = copy(dem)

    initial(a) = (true, typemax(a), a, limit)  # (is_pit, min, center, limit)
    @inbounds function update(v, a, _)
        if v[3] == a
            return v
        elseif (a - v[3]) > v[4]
            return (v[1] & true, min(a, v[2]), v[3], v[4])
        else
            (false, v[2], v[3], v[4])
        end
    end
    store!(d, i, v) = @inbounds d[i] = v[1] ? v[2] : v[3]

    return localfilter!(dst, dem, 3, initial, update, store!)
end

"""
    slope(dem::Matrix{<:Real}, cellsize=1.0, method=Horn())

Slope is the rate of change between a cell and its neighbors as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function slope(dem::AbstractMatrix{<:Real}; cellsize=1.0, method::DerivativeMethod=Horn())
    dst = copy(dem)
    slope(method, dst, dem, cellsize)
end


function slope(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] = atand(
        √(
            ((v[1] - v[2]) / (8 * v[5]))^2 + ((v[3] - v[4]) / (8 * v[5]))^2
        ))
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function slope(::ZevenbergenThorne, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] = atand(
        √(
            (v[1] / (2 * v[3]))^2 + (v[2] / (2 * v[3]))^2
        ))
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

function slope(::MDG, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        m = max(
            abs(v[1]) / (v[5] * sqrt2),
            abs(v[2]) / v[5],
            abs(v[3]) / (v[5] * sqrt2),
            abs(v[4]) / v[5],
        ) / 2
        d[i] = atand(√(m^2 + m^2))
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end


"""
    aspect(dem::Matrix{<:Real}, method=Horn())

Aspect is direction of [`slope`](@ref), as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function aspect(dem::AbstractMatrix{<:Real}; method::DerivativeMethod=Horn())
    dst = copy(dem)
    aspect(method, dst, dem, 1)
end

# Useless, as there's no x/y component.
function aspect(::MDG, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = max(
            abs(v[1]) / (v[5] * sqrt2),
            abs(v[2]) / v[5],
            abs(v[3]) / (v[5] * sqrt2),
            abs(v[4]) / v[5],
        ) / 2
        d[i] = compass(atand(δzδx, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end

function aspect(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = (v[1] - v[2]) / (8 * v[5])
        δzδy = (v[3] - v[4]) / (8 * v[5])
        # @info δzδx, δzδy
        d[i] = compass(atand(-δzδx, δzδy))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function aspect(::ZevenbergenThorne, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = v[1] / 2
        δzδy = v[2] / 2
        # @info δzδx, δzδy
        d[i] = compass(atand(δzδy, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

@inline @inbounds function zevenbergenthorne(v, a, b)
    if b == 2
        return (v[1], v[2] + a, v[3])
    elseif b == 4
        return (v[1] + a, v[2], v[3])
    elseif b == 6
        return (v[1] - a, v[2], v[3])
    elseif b == 8
        return (v[1], v[2] - a, v[3])
    else
        return v
    end
end

@inline @inbounds function horn(v, a, b)
    if b == 1
        return (v[1], v[2] + a, v[3], v[4] + a, v[5])
    elseif b == 2
        return (v[1], v[2] + 2a, v[3], v[4], v[5])
    elseif b == 3
        return (v[1], v[2] + a, v[3] + a, v[4], v[5])
    elseif b == 4
        return (v[1], v[2], v[3], v[4] + 2a, v[5])
    elseif b == 5
        return v
    elseif b == 6
        return (v[1], v[2], v[3] + 2a, v[4], v[5])
    elseif b == 7
        return (v[1] + a, v[2], v[3], v[4] + a, v[5])
    elseif b == 8
        return (v[1] + 2a, v[2], v[3], v[4], v[5])
    elseif b == 9
        return (v[1] + a, v[2], v[3] + a, v[4], v[5])
    end
end

@inbounds function mdg(v, a, b)
    if b == 1
        return (v[1] + a, v[2], v[3], v[4], v[5])
    elseif b == 2
        return (v[1], v[2], v[3], v[4] + a, v[5])
    elseif b == 3
        return (v[1], v[2], v[3] - a, v[4], v[5])
    elseif b == 4
        return (v[1], v[2] + a, v[3], v[4], v[5])
    elseif b == 5
        return v
    elseif b == 6
        return (v[1], v[2] - a, v[3], v[4], v[5])
    elseif b == 7
        return (v[1], v[2], v[3] + a, v[4], v[5])
    elseif b == 8
        return (v[1], v[2], v[3], v[4] - a, v[5])
    elseif b == 9
        return (v[1] - a, v[2], v[3], v[4], v[5])
    end
end

function compass(aspect)
    a = 90 - aspect
    a < 0 && (a += 360)
    return a
end

function aspect(compass::Real)
    return (450 - compass) % 360
end

"""
    curvature(dem::Matrix{<:Real})

Curvature is derivative of [`slope`](@ref), so the second derivative of the DEM.
"""
function curvature(dem::AbstractMatrix{<:Real}; cellsize=1.0)
    dst = copy(dem)

    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function curvature(v, a, b)
        if b == 2
            return (v[1], v[2] + a, v[3], v[4])
        elseif b == 4
            return (v[1] + a, v[2], v[3], v[4])
        elseif b == 5
            return (v[1], v[2], v[3] + a, v[4])
        elseif b == 6
            return (v[1] + a, v[2], v[3], v[4])
        elseif b == 8
            return (v[1], v[2] + a, v[3], v[4])
        else
            return v
        end
    end
    function store!(d, i, v)
        δzδx = (v[1] / 2 - v[3]) / v[4]^2
        δzδy = (v[2] / 2 - v[3]) / v[4]^2
        d[i] = -2(δzδx + δzδy) * 100
    end
    return localfilter!(dst, dem, nbkernel, initial, curvature, store!)
end


"""
    hillshade(dem::Matrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=1.0)

hillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref) given a light source with azimuth and zenith angles in °, , as defined in
Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function hillshade(dem::AbstractMatrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=1.0)
    dst = similar(dem, UInt8)
    zenithr = deg2rad(zenith)
    azimuthr = deg2rad(aspect(azimuth))

    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
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
        d[i] = round(UInt8, max(0, 255 * ((cos(zenithr) * cos(slope)) + (sin(zenithr) * sin(slope) * cos(azimuthr - a)))))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

"""
    multihillshade(dem::Matrix{<:Real}; cellsize=1.0)

multihillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref). Like [`hillshade`](@ref), but now using multiple sources as defined in
https://pubs.usgs.gov/of/1992/of92-422/of92-422.pdf, similar to GDALs -multidirectional.
"""
function multihillshade(dem::AbstractMatrix{<:Real}; cellsize=1.0)
    dst = similar(dem, UInt8)
    zenithr = deg2rad(60)

    initial(A) = (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
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
        something = (
            w225 * (α + β * cos(deg2rad(aspect(225)) - a)) +
            w270 * (α + β * cos(deg2rad(aspect(270)) - a)) +
            w315 * (α + β * cos(deg2rad(aspect(315)) - a)) +
            w360 * (α + β * cos(deg2rad(aspect(360)) - a))) / 2

        d[i] = round(UInt8, max(0, 255 * something))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end
