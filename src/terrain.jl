const neib_8 = @SMatrix [1.0 1 1; 1 0 1; 1 1 1]
const neib_8_inf = @SMatrix [1.0 1 1; 1 Inf 1; 1 1 1]
const neib_8_mask = @SMatrix Bool[1 1 1; 1 0 1; 1 1 1]

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
    oA[begin, begin] = A[begin, begin]
    oA[begin, end] = A[begin, end]
    oA[end, begin] = A[end, begin]
    oA[end, end] = A[end, end]
    return oA
end

function terrain_kernel(dem::AbstractMatrix{<:Real}, f::Function)
    dembuffer = buffer_array(dem)
    output = similar(dem)
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
    terrain_kernel(dem, roughness_kernel)
end

function roughness_kernel(A)
    x = A .- A[2, 2]
    return maximum(abs.(x))
end

"""
    TPI(dem::Matrix{<:Real})

TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).
"""
function TPI(dem::AbstractMatrix{<:Real})
    terrain_kernel(dem, TPI_kernel)
end

function TPI_kernel(A)
    x = A .* neib_8
    return A[2, 2] - sum(x) / 8
end

"""
    TRI(dem::Matrix{<:Real})

TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells.
This algorithm uses the square root of the sum of the square of the difference between a central pixel and its surrounding cells.
This is recommended for terrestrial use cases.
"""
function TRI(dem::AbstractMatrix{<:Real})
    terrain_kernel(dem, TRI_kernel)
end

function TRI_kernel(A)
    x = (A .- A[2, 2]) .^ 2
    return sqrt(sum(x))
end


"""
    pitremoval(dem::Matrix{<:Real})

Remove pits from a DEM Array if the center cell of a 3x3 patch is `limit` lower or than the surrounding cells.
"""
function pitremoval(dem::AbstractMatrix{<:Real}, limit=2.0)
    terrain_kernel(dem, f -> _pitremoval(f, limit))
end

function _pitremoval(A, limit)
    order = sortperm(vec(A))
    return (order[1] == 5) && ((A[order[2]] - A[order[1]]) >= limit) ? Inf : A[2, 2]
end

"""
    slope(dem::Matrix{<:Real})

Slope is the rate of change between a cell and its neighbors as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function slope(dem::AbstractMatrix{<:Real}, cellsize=1.0)
    terrain_kernel(dem, f -> slope_kernel(f, cellsize))
end

"""
    aspect(dem::Matrix{<:Real})

Aspect is direction of [`slope`](@ref), as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function aspect(dem::AbstractMatrix{<:Real})
    terrain_kernel(dem, aspect_kernel)
end

function slope_kernel(A, cellsize=1.0)
    δzδx, δzδy = derivative(A, cellsize)
    return atand(√(δzδx^2 + δzδy^2))
end

function aspect_kernel(A)
    δzδx, δzδy = derivative(A)
    return compass(atand(δzδy, -δzδx))
end

function derivative(A, cellsize=1.0)
    δzδx = ((A[1, 3] + 2A[2, 3] + A[3, 3]) - (A[1, 1] + 2A[2, 1] + A[3, 1])) / (8 * cellsize)
    δzδy = ((A[3, 1] + 2A[3, 2] + A[3, 3]) - (A[1, 1] + 2A[1, 2] + A[1, 3])) / (8 * cellsize)
    return δzδx, δzδy
end

function compass(aspect)
    a = 90 - aspect
    a < 0 && (a += 360)
    return a
end
