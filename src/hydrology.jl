
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

nbs =
    CartesianIndex.([(-1, -1), (-1, 1), (1, -1), (1, 1), (-1, 0), (0, -1), (0, 1), (1, 0)])

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

"""
    _orient(ci::CartesianIndex{2}, cellsize)

Convert a pixel CartesianIndex offset to the table convention (dim1=East+, dim2=North+).
Accounts for the sign of `cellsize`: a GeoTIFF typically has negative `cellsize[2]`
(+dim2 = South), so the second component is flipped. This function is its own inverse.
"""
@inline function _orient(ci::CartesianIndex{2}, cellsize)
    i, j = Tuple(ci)
    CartesianIndex(i * Int(sign(cellsize[1])), j * Int(sign(cellsize[2])))
end

# Neighbor offsets in table convention (dim1=East+, dim2=North+), ordered by compass bearing
nbb2 =
    CartesianIndex.([
        (0, 1),    # N  - 0°
        (1, 1),    # NE - 45°
        (1, 0),    # E  - 90°
        (1, -1),   # SE - 135°
        (0, -1),   # S  - 180°
        (-1, -1),  # SW - 225°
        (-1, 0),   # W  - 270°
        (-1, 1),   # NW - 315°
        (0, 1),    # N  - 0°   (wrap)
        (1, 1),    # NE - 45°  (wrap)
        (1, 0),    # E  - 90°  (wrap)
        (1, -1),   # SE - 135° (wrap)
    ])

function infc(aspect, cellsize)
    # Normalize compass bearing to 0-360 range (0° = North, clockwise)
    aspect = mod(aspect, 360)

    # Each direction covers a 45° sector
    sector = floor(Int, aspect / 45)

    # Get the two neighboring directions (in table convention)
    dir1_idx = sector + 1
    dir2_idx = (sector + 1) % 8 + 1

    # Convert from table to pixel convention
    return _orient(nbb2[dir1_idx], cellsize), _orient(nbb2[dir2_idx], cellsize)
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

"""
    flowaccumulation(dem::AbstractMatrix, closed::Matrix{Bool}, method::FlowDirectionMethod)

Computes the flow accumulation of a digital elevation model (DEM) `dem` with an optional `closed` mask and a `method` for flow direction.
Returns the flow accumulation and the flow direction (local drainage direction or ldd)
"""
function flowaccumulation(
    dem::AbstractMatrix,
    closed = falses(size(dem));
    method = DInf(),
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
    method = DInf(),
    cellsize = cellsize(dem),
)
    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    order = ones(Int64, length(closed) - sum(closed))

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

    dirs = _accumulate!(method, acc, order, dir, R, dem, cellsize)
    return acc, dirs
end

function _accumulate!(::D8, acc, order, dir, R, dem, cellsize)
    for i in reverse(order)
        dir[i] == CartesianIndex(0, 0) && continue
        acc[R[i] + dir[i]] += acc[i]
    end
    output = similar(dem, FlowDirection{LDD, UInt8})
    output .= getindex.(Ref(_ldd_ci2dir), _orient.(dir, Ref(cellsize)))
    return output
end
function _accumulate!(::DInf, acc, order, dir, R, dem, cellsize)
    asp = aspect(dem; method = Horn(), cellsize = abs.(cellsize))
    visited = falses(size(acc))
    output = similar(dem, FlowDirection{D8D, UInt8})
    fill!(output, 0)

    for i in reverse(order)
        dir[i] == CartesianIndex(0, 0) && continue
        aspect = asp[i]

        if !isfinite(aspect)
            acc[R[i] + dir[i]] += acc[i]
            output[i] = _d8_ci2dir[_orient(dir[i], cellsize)]
            visited[i] = true
            continue
        end

        a, b = infc(aspect, cellsize)
        aa, bb = infa(aspect)

        # Depression
        if (a != dir[i] && b != dir[i])
            acc[R[i] + dir[i]] += acc[i]
            output[i] = _d8_ci2dir[_orient(dir[i], cellsize)]
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

        dirs = zero(UInt8)
        if R[i] + a in R && aa > 0
            acc[R[i] + a] += acc[i] * aa
            dirs |= _d8_ci2dir[_orient(a, cellsize)]
        end
        if R[i] + b in R && bb > 0
            acc[R[i] + b] += acc[i] * bb
            dirs |= _d8_ci2dir[_orient(b, cellsize)]
        end
        output[i] = dirs
        if visited[R[i] + a] && visited[R[i] + b]
            error()
            acc[R[i] + dir[i]] += acc[i]
        end
        visited[i] = true
    end
    return output
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
    nb = vec(collect(CartesianIndices(dists)) .- CartesianIndex(2, 2))
    output = similar(dem, FlowDirection{D8D, UInt8})
    fill!(output, 0)

    weights = zeros(size(contour_lengths))
    Σw = 0.0
    for i in reverse(order)
        dir[i] == CartesianIndex(0, 0) && continue
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
        if iszero(Σw) || iszero(weights[CartesianIndex(2, 2) + dir[i]])
            acc[R[i] + dir[i]] += acc[i]
            output[i] = _d8_ci2dir[_orient(dir[i], cellsize)]
            visited[i] = true
            continue
        end
        dirs = zero(UInt8)
        for (ri, weight) in enumerate(weights)
            iszero(weight) && continue
            I = R[i] + nb[ri]
            I in R || continue
            acc[I] += acc[i] * (weight / Σw)
            dirs |= _d8_ci2dir[_orient(nb[ri], cellsize)]
        end
        output[i] = dirs
        visited[i] = true
        fill!(weights, 0)
    end
    return output
end

"""
    topographic_wetness_index(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))

Computes the Topographic Wetness Index (TWI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function topographic_wetness_index(
    dem::AbstractMatrix;
    method = DInf(),
    cellsize = cellsize(dem),
)
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc / tand(s))
end
@deprecate TWI topographic_wetness_index

"""
    stream_power_index(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))

Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function stream_power_index(dem::AbstractMatrix; method = DInf(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end
@deprecate SPI stream_power_index

"""
    height_above_nearest_drainage(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem), threshold=1e10)

Compute Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem` 
with an optional `method` for flow direction, a `cellsize`, and an flowaccumulation `threshold` for stream definition.
"""
function height_above_nearest_drainage(
    dem::AbstractMatrix;
    method = DInf(),
    cellsize = cellsize(dem),
    threshold = 100,
)
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

"""
    height_above_nearest_drainage(dem::AbstractMatrix, stream_mask::AbstractMatrix{Bool})

Computes the Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem` 
given a stream definition as a boolean `stream_mask`.
"""
function height_above_nearest_drainage(
    dem::AbstractMatrix,
    stream_mask::AbstractMatrix{Bool},
)
    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    closed = falses(size(dem))
    order = ones(Int64, length(closed) - sum(closed))

    output = zero(dem)

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

"""
    depression_depth(dem::AbstractMatrix; filled=filldepressions(dem))

Computes the depth of each cell below the filled surface.

Returns the difference between the depression-filled DEM and the original DEM,
representing how deep each cell sits within a depression/depression. Cells not in
depressions will have a depth of zero.

This is useful for identifying potential cold air pooling zones and water
retention areas.
"""
function depression_depth(dem::AbstractMatrix; filled = filldepressions(dem))
    filled .- dem
end

"""
    depression_volume(dem::AbstractMatrix; filled=filldepressions(dem), cellsize=cellsize(dem))

Computes the total volume of all depressions/basins in the DEM.

Returns the sum of depression depths multiplied by cell area.
"""
function depression_volume(
    dem::AbstractMatrix;
    filled = filldepressions(dem),
    cellsize = cellsize(dem),
)
    sum(filter(isfinite, depression_depth(dem; filled))) * prod(abs.(cellsize))
end

"""
    drainage_potential(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))

Computes a drainage potential index indicating how well each cell drains.

High values indicate good drainage (steep slopes with low upstream accumulation).
Low values indicate poor drainage (flat areas that accumulate flow from upslope).

This is useful as a proxy for cold air drainage - cells with high drainage potential
will shed cold air downslope rather than pooling it.

Based on the relationship between slope and flow accumulation:
`drainage = sin(slope) / log(1 + accumulation)`
"""
function drainage_potential(dem::AbstractMatrix; method = DInf(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. sind(s) / log1p(acc)
end
