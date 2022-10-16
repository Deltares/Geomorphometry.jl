const neib_8 = @SMatrix [1.0 1 1; 1 0 1; 1 1 1]
const neib_8_inf = @SMatrix [1.0 1 1; 1 Inf 1; 1 1 1]
const neib_8_mask = @SMatrix Bool[1 1 1; 1 0 1; 1 1 1]

abstract type DerivativeMethod end

"""Second order finite difference estimator using all 4 neighbors (Zevenbergen and Thorne, 1987)."""
struct ZevenbergenThorne <: DerivativeMethod end

"""Third order finite difference estimator using all 8 neighbors (Horn, 1981)."""
struct Horn <: DerivativeMethod end

"""Maximum Downward Gradient"""
struct MDG <: DerivativeMethod end

function noise(factor=0.01)
    (rand() - 0.5) * factor
end

function buffer_array(A::AbstractMatrix{<:Real})
    oA = OffsetMatrix(fill(zero(eltype(A)), size(A) .+ 2), UnitRange.(0, size(A) .+ 1))
    # Update center
    oA[begin+1:end-1, begin+1:end-1] .= A
    # Set edges to mirror center
    oA[begin, begin+1:end-1] .= A[begin, :]
    oA[end, begin+1:end-1] .= A[end, :]
    oA[begin+1:end-1, begin] .= A[:, begin]
    oA[begin+1:end-1, end] .= A[:, end]
    # Set corners to mirror corners of center
    oA[begin, begin] = A[begin, begin] + noise()
    oA[begin, end] = A[begin, end] + noise()
    oA[end, begin] = A[end, begin] + noise()
    oA[end, end] = A[end, end] + noise()
    return oA
end

function terrain_kernel(dem::AbstractMatrix{T}, f::Function, t=T) where {T<:Real}
    dembuffer = buffer_array(dem)
    output = similar(dem, t)
    Δ = CartesianIndex(1, 1)

    @inbounds for I ∈ CartesianIndices(size(output))
        patch = I-Δ:I+Δ
        dempatch = SMatrix{3,3}(view(dembuffer, patch))
        output[I] = f(dempatch)
    end
    output
end

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
    terrain_kernel(dem, f -> slope_kernel(f, cellsize, method))
end

function slope_kernel(A, cellsize=1.0, method::DerivativeMethod=Horn())
    δzδx, δzδy = derivative(method, A, cellsize)
    return atand(√(δzδx^2 + δzδy^2))
end

"""
    aspect(dem::Matrix{<:Real}, method=Horn())

Aspect is direction of [`slope`](@ref), as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function aspect(dem::AbstractMatrix{<:Real}; method::DerivativeMethod=Horn())
    terrain_kernel(dem, f -> aspect_kernel(f, method))
end

function aspect_kernel(A, method::DerivativeMethod=Horn())
    δzδx, δzδy = derivative(method, A, 1.0)
    return compass(atand(δzδy, -δzδx))
end

function derivative(::ZevenbergenThorne, A, cellsize=1.0)
    δzδx = (A[1, 2] - A[3, 2]) / 2 * cellsize
    δzδy = (A[2, 1] - A[2, 3]) / 2 * cellsize
    return δzδx, δzδy
end

function derivative(::Horn, A, cellsize=1.0)
    δzδx = ((A[1, 3] + 2A[2, 3] + A[3, 3]) - (A[1, 1] + 2A[2, 1] + A[3, 1])) / (8 * cellsize)
    δzδy = ((A[3, 1] + 2A[3, 2] + A[3, 3]) - (A[1, 1] + 2A[1, 2] + A[1, 3])) / (8 * cellsize)
    return δzδx, δzδy
end

function derivative(::MDG, A, cellsize=1.0)
    δzδx = max(
        abs(A[1, 1] - A[3, 3]) / (cellsize * sqrt2),
        abs(A[1, 2] - A[3, 2]) / cellsize,
        abs(A[1, 3] - A[3, 1]) / (cellsize * sqrt2),
        abs(A[2, 1] - A[2, 3]) / cellsize,
    ) / 2

    return δzδx, δzδx
end

function compass(aspect)
    a = 90 - aspect
    a < 0 && (a += 360)
    return a
end

function aspect(compass::Real)
    return (360 - compass + 90) % 360
end

"""
    curvature(dem::Matrix{<:Real})

Curvature is derivative of [`slope`](@ref), so the second derivative of the DEM.
"""
function curvature(dem::AbstractMatrix{<:Real}; cellsize=1.0)
    return terrain_kernel(dem, f -> curvature_kernel(f, cellsize))
end

function curvature_kernel(A, cellsize=1.0)
    δzδx = ((A[1, 2] + A[3, 2]) / 2 - A[2, 2]) / cellsize^2
    δzδy = ((A[2, 1] + A[2, 3]) / 2 - A[2, 2]) / cellsize^2

    return -2(δzδx + δzδy) * 100
end

"""
    hillshade(dem::Matrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=1.0, method=Horn())

hillshade is the simulated illumination of a surface based on its [`slope`](@ref) and
[`aspect`](@ref) given a light source with azimuth and zenith angles in °, , as defined in
Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function hillshade(dem::AbstractMatrix{<:Real}; azimuth=315.0, zenith=45.0, cellsize=1.0, method::DerivativeMethod=Horn())
    return terrain_kernel(dem, f -> hillshade_kernel(f, cellsize, azimuth, zenith, method), UInt8)
end

function hillshade_kernel(A, cellsize=1.0, azimuth=315.0, zenith=45.0, method::DerivativeMethod=Horn())
    z = deg2rad(zenith)
    c = deg2rad(aspect(azimuth))

    δzδx, δzδy = derivative(method, A, cellsize)
    if δzδx != 0
        a = atan(δzδy, -δzδx)
        if a < 0
            a += 2π
        end
    else
        a = π / 2
        if δzδy > 0
            a += 2π
        end
    end

    s = atan(√(δzδx^2 + δzδy^2))
    return round(UInt8, 255 * ((cos(z) * cos(s)) + (sin(z) * sin(s) * cos(c - a))))
end
