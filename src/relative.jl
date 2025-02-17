"""
    roughness(dem::AbstractMatrix{<:Real})

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
    TPI(dem::AbstractMatrix{<:Real})

TPI stands for Topographic Position Index, which is defined as the difference between a central pixel and the mean of its surrounding cells (see Wilson et al 2007, Marine Geodesy 30:3-35).
"""
function TPI(dem::AbstractMatrix{<:Real}, window::Stencil = Moore(1))
    mapstencil(x -> center(x) - mean(x), window, dem)
end

"""
    BPI(dem::AbstractMatrix{<:Real})

BPI stands for Bathymetric Position Index (Lundblad et al., 2006), which is defined as the difference between a central pixel and the mean of the cells in an annulus around it.
"""
function BPI(dem::AbstractMatrix{<:Real}, window::Annulus = Annulus(3, 2))
    mapstencil(x -> center(x) - mean(x), window, dem)
end

"""
    RIE(dem::AbstractMatrix{<:Real})

RIE stands for Roughness Index Elevation, which quantifies the standard deviation of residual topography (Cavalli et al., 2008)
"""
function RIE(dem::AbstractMatrix{<:Real}, window::Stencil = Window(1))
    meandem = mapstencil(x -> center(x) - mean(x), window, dem)
    mapstencil(std, window, meandem)
end

"""
    TRI(dem::AbstractMatrix{<:Real})

TRI stands for Terrain Ruggedness Index, which measures the difference between a central pixel and its surrounding cells.
This algorithm uses the square root of the sum of the square of the difference between a central pixel and its surrounding cells.
This is recommended for terrestrial use cases.
"""
function TRI(dem::AbstractMatrix{<:Real}; normalize = false, squared = true)
    dst = copy(dem)

    @inline initial(a) = (zero(a), a)
    @inline update(v, a, _) =
        if squared
            (v[1] + (a - v[2])^2, v[2])
        else
            (v[1] + abs(a - v[2]), v[2])
        end
    @inline store!(d, i, v) =
        if normalize && squared
            @inbounds d[i] = sqrt(v[1] / 8)
        elseif normalize
            @inbounds d[i] = v[1] / 8
        elseif squared
            @inbounds d[i] = sqrt(v[1])
        else
            @warn "TRI: normalize=false and squared=false is not recommended."
            @inbounds d[i] = v[1]
        end

    return localfilter!(dst, dem, 3, initial, update, store!)
end

"""
    prominence(dem::AbstractMatrix{<:Real})

Prominence calculates the number of cells that are lower or equal than the central cell.
Thus, 8 is a local maximum (peak), while 0 is a local minimum (pit).
"""
function prominence(dem::AbstractMatrix{<:Real})
    dst = similar(dem, Int8)

    initial(a) = (; count = 0, center = a)
    update(v, a, _) = (; count = v.count + (a <= v.center), center = v.center)
    store!(d, i, v) = @inbounds d[i] = v.count - 1

    return localfilter!(dst, dem, 3, initial, update, store!)
end

round_step(x, step) = round(x / step) * step

"""
    entropy(dem::AbstractMatrix{<:Real}; step=0.5)

Entropy calculates the Shannon entropy of the surrounding cells of a central cell.
`step` is the bin size for the histogram.
"""
function entropy(dem::AbstractMatrix{<:Real}; step = 0.5)
    if !isnothing(step)
        dem = round_step.(dem, step)
    end
    dst = copy(dem)
    entropy!(dst, dem)
end

"""
    entropy!(dem::AbstractMatrix{<:Real})

Entropy calculates the Shannon entropy of the surrounding cells of a central cell.
"""
function entropy!(dst::AbstractMatrix{<:Real}, dem::AbstractMatrix{T}) where {T <: Real}

    # TODO Exclude center cell?

    # Manually setup buffers to avoid allocations.
    values = @MVector zeros(T, 9)
    counts = @MVector zeros(Float32, 9)  # Float32 so it can be normalized in place.

    # As previously set values are not removed, start at 0 
    # to always have a new value at the first iteration.
    initial(_) = (; counts, values, i = 0)
    function update(v, a, _)
        # We either match an existing value
        for i in 1:(v.i)
            if a == v.values[i]
                v.counts[i] += 1
                return v
            end
        end
        # Or add a new value
        v.values[v.i + 1] = a
        v.counts[v.i + 1] = 1
        return (; counts = v.counts, values = v.values, i = v.i + 1)
    end
    function store!(d, i, v)
        @inbounds d[i] = @views _entropy(v.counts[begin:(v.i)], v.values[begin:(v.i)])
        v.counts .= 0  # reset counts for next iteration
    end

    return localfilter!(dst, dem, 5, initial, update, store!)
end

function _entropy(counts, values)
    # Counts and values are re-used for the probability and log(probability) calculations.
    prob = counts ./= sum(counts)
    values .= log.(prob)
    prob .*= values
    return -sum(prob)
end

function cross2(a, b)
    if !(length(a) == length(b) == 3)
        throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
    end
    a1, a2, a3 = a
    b1, b2, b3 = b
    return (a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)
end

"""
    rugosity(dem::AbstractMatrix{<:Real})

Compute the rugosity of a DEM, which is the ratio between the 
surface area divided by the planimetric area.

Jenness 2019 https://onlinelibrary.wiley.com/doi/abs/10.2193/0091-7648%282004%29032%5B0829%3ACLSAFD%5D2.0.CO%3B2
"""
function rugosity(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem))
    dst = zeros(Float32, size(dem))

    δx, δy = cellsize

    # Manually setup buffers to avoid allocations.
    values = @MVector zeros(Float32, 9)
    mask = falses(9)

    initial(a) = (; values = values, mask = mask, initial = a)
    function update(v, a, b)
        v.values[b] = a - v.initial
        v.mask[b] = true
        return (; values = v.values, mask = mask, initial = v.initial)
    end
    function store!(d, i, v)
        d[i] += area((δx, δy, v.values[1]), (δx, 0.0, v.values[2]))
        d[i] += area((δx, δy, v.values[3]), (δx, 0.0, v.values[2]))
        d[i] += area((δx, δy, v.values[3]), (δy, 0.0, v.values[6]))
        d[i] += area((δx, δy, v.values[9]), (δy, 0.0, v.values[6]))
        d[i] += area((δx, δy, v.values[9]), (δx, 0.0, v.values[8]))
        d[i] += area((δx, δy, v.values[7]), (δx, 0.0, v.values[8]))
        d[i] += area((δx, δy, v.values[7]), (δy, 0.0, v.values[4]))
        d[i] += area((δx, δy, v.values[1]), (δy, 0.0, v.values[4]))
        fill!(v.values, 0)
        fill!(v.mask, false)
        d[i] /= (δx * δy)
    end

    return localfilter!(dst, dem, nbkernel, initial, update, store!)
end

function area(a, b)
    # Divide by 2 for the triangle area, divide by 4 to get the area within the centre cell
    sqrt(sum(cross2(a, b) .^ 2)) / 8
end

"""
    pitremoval(dem::AbstractMatrix{<:Real})

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
