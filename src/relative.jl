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
function TPI(dem::AbstractMatrix{<:Real}, window::Stencil = Moore(1))
    mapstencil(x -> center(x) - mean(x), window, dem)
end

"""
    BPI(dem::Matrix{<:Real})

BPI stands for Bathymetric Position Index (Lundblad et al., 2006), which is defined as the difference between a central pixel and the mean of the cells in an annulus around it.
"""
function BPI(dem::AbstractMatrix{<:Real}, window::Annulus = Annulus(3, 2))
    mapstencil(x -> center(x) - mean(x), window, dem)
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
    RIE(dem::Matrix{<:Real})

RIE stands for Roughness Index Elevation, which quantifies the standard deviation of residual topography (Cavalli et al., 2008)
"""
function RIE(dem::AbstractMatrix{<:Real}, window::Stencil = Window(1))
    meandem = mapstencil(x -> center(x) - mean(x), window, dem)
    mapstencil(std, window, meandem)
end

"""
    pitremoval(dem::Matrix{<:Real})

Remove pits from a DEM Array if the center cell of a 3x3 patch is `limit` lower or than the surrounding cells.
"""
function pitremoval(dem::AbstractMatrix{<:Real}; limit = 0.0)
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
