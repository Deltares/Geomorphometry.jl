
function edges(A::AbstractMatrix)
    CI = CartesianIndices(A)
    edges(CI)
end

function edges(CI::CartesianIndices)
    indices = Vector{CartesianIndex}()
    append!(indices, first(eachrow(CI)))
    append!(indices, last(eachrow(CI)))
    append!(indices, first(eachcol(CI)))
    append!(indices, last(eachcol(CI)))
    unique(indices)
end

"""
    filldepressions(dem::AbstractMatrix, mask::Matrix{Bool})

Performs the Priority-Flood algorithm [barnesPriorityFloodOptimalDepressionFilling2014](@cite) on the given digital elevation model (DEM) `dem` with an optional `mask`.

# Arguments
- `dem::AbstractMatrix`: A 2D array representing the digital elevation model (DEM).
- `mask::AbstractMatrix{Bool}`: A 2D boolean array representing the mask. Cells with `true` values are considered in the computation, while cells with `false` values are ignored.
"""
function filldepressions(dem::AbstractMatrix, mask = falses(size(dem)))
    filldepressions!(copy(dem), mask)
end

const Δ = CartesianIndex(1, 1)

abstract type FlowDirectionMethod end

"""D8 Flow Direction method by [Jenson (1988)](@cite jensonExtractingTopographicStructure1988)."""
struct D8 <: FlowDirectionMethod end

"""DInf Flow Direction method by [Tarboton (1997)](@cite tarbotonNewMethodDetermination1997)."""
struct DInf <: FlowDirectionMethod end

"""FD8 Flow Direction method by [Quin (1991)](@cite quinnPredictionHillslopeFlow1991)."""
Base.@kwdef struct FD8 <: FlowDirectionMethod
    p::Float32 = 1.1
end

function filldepressions!(dem::AbstractMatrix, queued = falses(size(dem)))
    open = PriorityQueue{CartesianIndex{2}, eltype(dem)}()
    pit = DataStructures.Queue{CartesianIndex{2}}()

    R = CartesianIndices(dem)
    I_first, I_last = first(R), last(R)

    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        queued[cell] = true  # queued
    end
    @inbounds while !isempty(open) || !isempty(pit)
        cell = !isempty(pit) ? DataStructures.dequeue!(pit) : dequeue!(open)
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true
            if dem[ncell] <= dem[cell]
                dem[ncell] = dem[cell]
                DataStructures.enqueue!(pit, ncell)
            else
                enqueue!(open, ncell, dem[ncell])
            end
        end
    end
    return dem
end

nbs = CartesianIndex.([
    (-1, -1),
    (-1, 1),
    (1, -1),
    (1, 1),
    (-1, 0),
    (0, -1),
    (0, 1),
    (1, 0),
])

function watersheds(dem::AbstractMatrix, queued = falses(size(dem)))
    open = PriorityQueue{CartesianIndex{2}, eltype(dem)}()
    pit = DataStructures.Queue{CartesianIndex{2}}()
    labels = zeros(Int, size(dem))
    label = 1
    pits = falses(size(dem))

    R = CartesianIndices(dem)
    I_first, I_last = first(R), last(R)

    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        queued[cell] = true  # queued
    end
    @inbounds while !isempty(open) || !isempty(pit)
        cell = !isempty(pit) ? DataStructures.dequeue!(pit) : dequeue!(open)
        if queued[cell] &&
           !ismissing(dem[cell]) &&
           isfinite(dem[cell]) &&
           iszero(labels[cell])
            labels[cell] = label
            label += 1
        end
        # for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
        for nb in nbs
            ncell = cell + nb
            ncell in R || continue
            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true
            if dem[ncell] <= dem[cell]
                pits[ncell] = true
                labels[ncell] = label
                label += 1
                dem[ncell] = dem[cell]
                DataStructures.enqueue!(pit, ncell)
            else
                if !ismissing(dem[ncell]) && isfinite(dem[ncell])
                    if pits[cell]
                        labels[ncell] = label
                        label += 1
                    else
                        labels[ncell] = labels[cell]
                    end
                end
                enqueue!(open, ncell, dem[ncell])
            end
        end
    end
    return dem, labels
end

# 1  2  3
# 4  5  6
# 7  8  9

nbb2 =
    CartesianIndex.([
        (-1, 0),  # N  - 270°
        (-1, 1),  # NE - 315°
        (0, 1),   # E  - 0° 
        (1, 1),   # SE - 45°
        (1, 0),   # S  - 90°
        (1, -1),  # SW - 135°
        (0, -1),  # W  - 180°
        (-1, -1), # NW - 225°
        (-1, 0),  # N  - 270°
        (-1, 1),  # NE - 315°
        (0, 1),   # E  - 0° 
        (1, 1),   # SE - 45°
    ])

function infc(aspect)
    # Normalize aspect to 0-360 range
    aspect = mod(aspect+90, 360)
    
    # Convert aspect to index (0° = North, clockwise)
    # Each direction covers 45° sector
    sector = floor(Int, aspect / 45)
    
    # Get the two neighboring directions
    dir1_idx = sector + 1
    dir2_idx = (sector + 1) % 8 + 1
    
    return nbb2[dir1_idx], nbb2[dir2_idx]
end

function infa(aspect)
    # Normalize aspect to 0-360 range
    aspect = mod(aspect, 360)
    
    # Find position within 45° sector
    sector_angle = mod(aspect, 45)
    
    # Calculate weights for the two directions
    # Weight decreases linearly from 1 to 0 as we move away from the direction
    weight2 = sector_angle / 45
    weight1 = 1 - weight2
    
    return weight1, weight2
end

const directions = centered([1 2 3; 4 5 6; 7 8 9])

"""
    flowaccumulation(dem::AbstractMatrix, closed::Matrix{Bool}, method::FlowDirectionMethod)

Computes the flow accumulation of a digital elevation model (DEM) `dem` with an optional `closed` mask and a `method` for flow direction.
Returns the flow accumulation and the flow direction (local drainage direction or ldd)
"""
function flowaccumulation(
    dem::AbstractMatrix,
    closed = falses(size(dem));
    method = D8(),
    cellsize = cellsize(dem),
)
    acc = similar(dem, Float32)
    acc .= abs(cellsize[1] * cellsize[2])
    flowaccumulation!(dem, acc, copy(closed); method, cellsize)
end

function flowaccumulation!(
    dem::AbstractMatrix,
    acc::AbstractMatrix{<:Real},
    closed = falses(size(dem));
    method = D8(),
    cellsize = cellsize(dem),
)
    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    order = ones(Int64, length(closed) - sum(closed))

    output = similar(dem, eltype(directions))

    open = PriorityQueue{CartesianIndex{2}, eltype(dem)}()

    R = CartesianIndices(dem)
    L = LinearIndices(dem)
    I_first, I_last = first(R), last(R)

    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        # for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
        for nb in nbs
            ncell = cell + nb
            ncell in R || continue
            # skip visited and center cells
            (closed[ncell] || ncell == cell) && continue

            closed[ncell] = true
            dir[ncell] = cell - ncell

            enqueue!(open, ncell, dem[ncell])
        end
    end

    _accumulate!(method, acc, order, dir, R, dem, cellsize)
    output .= getindex.(Ref(directions), dir)
    return acc, output
end

function _accumulate!(::D8, acc, order, dir, R, dem, cellsize)
    for i in reverse(order)
        acc[R[i] + dir[i]] += acc[i]
    end
end
function _accumulate!(::DInf, acc, order, dir, R, dem, cellsize)
    asp = aspect(dem; method = Horn(), cellsize=abs.(cellsize))
    visited = falses(size(acc))

    for i in reverse(order)
        aspect = asp[i]

        if !isfinite(aspect)
            acc[R[i] + dir[i]] += acc[i]
            visited[i] = true
            continue
        end

        a, b = infc(aspect)
        aa, bb = infa(aspect)

        # a, b = lookup[dir[i] + CartesianIndex(2, 2)]

        # Depression
        if (a != dir[i] && b != dir[i])
            acc[R[i] + dir[i]] += acc[i]
            # acc[R[i] + dir[i]] = 5000
            visited[i] = true
            continue
        end

        # Scale flows correctly at the edges
        if !(R[i] + a in R) || visited[R[i] + a]
            aa = 0
            bb = 1
        end
        if !(R[i] + b in R) || visited[R[i] + b]
            aa = 1
            bb = 0
        end

        if R[i] + a in R
            acc[R[i] + a] += acc[i] * aa
        end
        if R[i] + b in R
            acc[R[i] + b] += acc[i] * bb
        end
        if visited[R[i] + a] && visited[R[i] + b]
            error()
            acc[R[i] + dir[i]] += acc[i]
        end
        visited[i] = true
    end
end

function _accumulate!(fd8::FD8, acc, order, dir, R, dem, cellsize)
    # Derive contour lengths, which is used to calculate the weights
    # Uses the algorithm by Quinn et al. (1991), L1=0.5 L2=0.354 for δx=δy=1
    # TODO Check whether just using the angles is enough.
    δx, δy = abs.(cellsize)
    δxy = sqrt((δx / 4)^2 + (δy / 4)^2)
    contour_lengths = @SMatrix [
        δxy δx δxy
        δy 0 δy
        δxy δx δxy
    ]

    δxy = sqrt(δx^2 + δy^2)
    dists = @SMatrix [
        δxy δx δxy
        δy 0 δy
        δxy δx δxy
    ]

    visited = falses(size(acc))
    nb = vec(collect(CartesianIndices(dists)).-CartesianIndex(2,2))

    weights = zeros(size(contour_lengths))
    Σw = 0.0
    for i in reverse(order)
        Σw = 0.0
        # TODO Fix this distance with actual distances
        for (ri, dist) in enumerate(dists)
            ri == 5 && continue
            I = R[i] + nb[ri]
            I in R || continue
            visited[I] && continue
            diff = dem[R[i]] - dem[I]
            if diff < 0 || isnan(diff) # neighbor is higher
                continue
            end
            # TODO Check whether this diff/dist is good enough
            weight = (diff / dist * contour_lengths[ri])^fd8.p
            weights[ri] = weight
            Σw += weight
        end
        if iszero(Σw) || iszero(weights[CartesianIndex(2,2)+dir[i]])
            acc[R[i] + dir[i]] += acc[i]
            visited[i] = true
            continue
        end
        for (ri, weight) in enumerate(weights)
            iszero(weight) && continue
            I = R[i] + nb[ri]
            I in R || continue
            acc[I] += acc[i] * (weight / Σw)
        end
        visited[i] = true
        fill!(weights, 0)
    end
end

"""
    TWI(dem::AbstractMatrix; method=D8(), cellsize=cellsize(dem))

Computes the Topographic Wetness Index (TWI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function TWI(dem::AbstractMatrix; method = D8(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc / tand(s))
end

"""
    SPI(dem::AbstractMatrix; method=D8(), cellsize=cellsize(dem))

Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function SPI(dem::AbstractMatrix; method = D8(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end

"""
    HAND(dem::AbstractMatrix; method=D8(), cellsize=cellsize(dem), threshold=1e10)

Computes the Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem` 
with an optional `method` for flow direction, a `cellsize`, and a `threshold` for stream definition.
"""
function HAND(dem::AbstractMatrix; method = D8(), cellsize = cellsize(dem), threshold = 100)
    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    closed = falses(size(dem))
    order = ones(Int64, length(closed) - sum(closed))

    output = zero(dem)
    acc = similar(dem, Float32)
    acc .= abs(cellsize[1] * cellsize[2])

    open = PriorityQueue{CartesianIndex{2}, eltype(dem)}()

    R = CartesianIndices(dem)
    L = LinearIndices(dem)
    I_first, I_last = first(R), last(R)

    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        # for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
        for nb in nbs
            ncell = cell + nb
            ncell in R || continue
            # skip visited and center cells
            (closed[ncell] || ncell == cell) && continue

            closed[ncell] = true
            dir[ncell] = cell - ncell

            enqueue!(open, ncell, dem[ncell])
        end
    end

    _accumulate!(method, acc, order, dir, R, dem, cellsize)
    stream_mask = acc .>= threshold
    _hand!(output, order, dir, R, dem, stream_mask)
    return output
end

function _hand!(output, order, dir, R, dem, stream_mask)
    for i in order
        if stream_mask[i]
            # Relative height for stream is 0
            output[i] = 0.0
        elseif isfinite(dem[i]) && isfinite(dem[R[i] + dir[i]])
            # Otherwise, add the height difference with the downstream cell
            # to the downstream cell's HAND value. Use max to avoid negative values
            # in case of depressions. TODO Handle depressions better?
            output[i] = output[R[i] + dir[i]] + max(0, (dem[i] - dem[R[i] + dir[i]]))
        else
            output[i] = NaN
        end
    end
    return output
end
