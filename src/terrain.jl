const neib_8 = @SMatrix[1. 1 1; 1 0 1; 1 1 1]

function buffer_array(A::Array)
    oA = OffsetMatrix(fill(Inf, size(A) .+ 2), UnitRange.(0, size(A) .+ 1))
    # Update center
    oA[begin + 1:end - 1,begin + 1:end - 1] .= A
    # Set edges to mirror center
    oA[begin, begin + 1:end - 1] .= A[begin, :]
    oA[end, begin + 1:end - 1] .= A[end, :]
    oA[begin + 1:end - 1, begin] .= A[:, begin]
    oA[begin + 1:end - 1, end] .= A[:, end]
    # Set corners to mirror corners of center
    oA[begin, begin] = A[begin, begin]
    oA[begin, end] = A[begin, end]
    oA[end, begin] = A[end, begin]
    oA[end, end] = A[end, end]
    return oA
end

"""
roughness(dem::Matrix{<:Real})

Roughness is the largest inter-cell difference of a central pixel and its surrounding cell, as defined in Wilson et al (2007, Marine Geodesy 30:3-35).
"""
function roughness(dem::Matrix{<:Real})

    ex_dem = buffer_array(dem)
    roughness = similar(dem)
    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(roughness))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= rdata .- rdata[2,2]
        roughness[I] = maximum(abs.(x))
    end
    roughness
end

"""
TPI(dem::Matrix{<:Real})

TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).
"""
function TPI(dem::Matrix{<:Real})

    ex_dem = buffer_array(dem)
    tpi = similar(dem)
    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(tpi))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= rdata .* neib_8
        tpi[I] = rdata[2,2] - mean(x)
    end
    return tpi
end

"""
TRI(dem::Matrix{<:Real})

TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells.
This algorithm uses the square root of the sum of the square of the difference between a central pixel and its surrounding cells.
This is recommended for terrestrial use cases.
"""
function TRI(dem::Matrix{<:Real})

    ex_dem = buffer_array(dem)
    tri = similar(dem)
    x = @MMatrix zeros(3, 3)

    @inbounds for I ∈ CartesianIndices(size(tri))
        patch = I - Δ:I + Δ
        rdata = view(ex_dem, patch)

        x .= (rdata .- rdata[2,2]).^2
        tri[I] = sqrt(sum(x))
    end
    return tri
end