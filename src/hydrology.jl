
function outlets(A::AbstractMatrix)
    CI = CartesianIndices(A)
    indices = Vector{CartesianIndex}()
    append!(indices, first(eachrow(CI)))
    append!(indices, last(eachrow(CI)))
    append!(indices, first(eachcol(CI)))
    append!(indices, last(eachcol(CI)))
    unique(indices)
end

_cells(dem) = eachindex(dem)
_cells(dem::AbstractMatrix) = CartesianIndices(dem)

cellarea(dem::AbstractMatrix, cell; cellsize = cellsize(dem)) =
    abs(cellsize[1] * cellsize[2])

"""
    filldepressions(dem, mask=fill!(similar(dem, Bool), false))

Performs the Priority-Flood algorithm [barnesPriorityFloodOptimalDepressionFilling2014](@cite) on the given digital elevation model (DEM) `dem` with an optional `mask`.

# Arguments
- `dem`: An indexed digital elevation model.
- `mask`: An optional boolean array representing the mask. Cells with `true` values are treated as already filled (queued) and excluded from the computation. Defaults to all `false`.
"""
function filldepressions(dem, mask = fill!(similar(dem, Bool), false))
    filldepressions!(copy(dem), mask)
end

abstract type FlowDirectionMethod end

"""D8 Flow Direction method by [Jenson (1988)](@cite jensonExtractingTopographicStructure1988)."""
struct D8 <: FlowDirectionMethod end

"""DInf Flow Direction method by [Tarboton (1997)](@cite tarbotonNewMethodDetermination1997)."""
struct DInf <: FlowDirectionMethod end

"""
    FD8(; p=1.1)

FD8 (multiple) Flow Direction method by [Quinn (1991)](@cite quinnPredictionHillslopeFlow1991).
The exponent `p` controls how strongly flow concentrates towards the steepest descent: higher
values approach single-direction flow, lower values spread flow more evenly.
"""
Base.@kwdef struct FD8 <: FlowDirectionMethod
    p::Float32 = 1.1
end

function filldepressions!(dem, queued = fill!(similar(dem, Bool), false))
    R = _cells(dem)
    first_cell = first(R)
    open = PriorityQueue{typeof(first_cell), eltype(dem)}()
    pit = DataStructures.Queue{typeof(first_cell)}()

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        queued[cell] = true  # queued
    end
    @inbounds while !isempty(open) || !isempty(pit)
        cell = !isempty(pit) ? DataStructures.dequeue!(pit) : dequeue!(open)
        for ncell in neighbors(dem, cell)
            ncell in R || continue
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

const nbs =
    CartesianIndex.(((-1, -1), (-1, 1), (1, -1), (1, 1), (-1, 0), (0, -1), (0, 1), (1, 0)))

"""
    neighbors(dem, cell)

Iterate over the neighbors of `cell` that fall within `dem`. Matrices use the
Moore neighborhood. Implementations must yield only in-domain indices: kernels
index `dem` with the result directly, some under `@inbounds`.
"""
neighbors(dem, cell::CartesianIndex{2}) =
    Iterators.filter(in(CartesianIndices(dem)), cell + nb for nb in nbs)

function watersheds(dem, queued = fill!(similar(dem, Bool), false))
    R = _cells(dem)
    first_cell = first(R)
    open = PriorityQueue{typeof(first_cell), eltype(dem)}()
    pit = DataStructures.Queue{typeof(first_cell)}()
    labels = zeros(Int, size(dem))
    label = 1
    pits = fill!(similar(dem, Bool), false)

    @inbounds for cell in outlets(dem)
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
        for ncell in neighbors(dem, cell)
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
@inline _orient(direction, cellsize) = direction

_relativeindices(::Type{R}, direction::FlowDirection, center, cellsize) where {R} =
    decompose(R, direction, center)
_relativeindices(
    ::Type{CartesianIndex{2}},
    direction::FlowDirection,
    center,
    cellsize,
) =
    (
        _orient(relative, cellsize) for
        relative in decompose(CartesianIndex{2}, direction)
    )

_directioncode(dem::AbstractMatrix, direction, cellsize) =
    _d8_ci2dir[_orient(direction, cellsize)]
function _directioncode(dem, direction, cellsize)
    UInt8(Int(FlowDirection{D8D}(direction)))
end

function _flowneighbors(dem, cell, aspect, cellsize)
    aspect = mod(aspect, 360)
    lower_cell = upper_cell = cell
    lower_gap = upper_gap = Inf
    found = false
    for neighbor in neighbors(dem, cell)
        found = true
        bearing = _cellbearing(dem, cell, neighbor, cellsize)
        candidate_lower_gap = mod(aspect - bearing, 360)
        candidate_upper_gap = mod(bearing - aspect, 360)
        if candidate_lower_gap < lower_gap
            lower_cell = neighbor
            lower_gap = candidate_lower_gap
        end
        if candidate_upper_gap < upper_gap
            upper_cell = neighbor
            upper_gap = candidate_upper_gap
        end
    end

    found || return nothing
    lower_cell == upper_cell && return lower_cell, upper_cell, 1.0, 0.0
    gap = lower_gap + upper_gap
    return lower_cell, upper_cell, upper_gap / gap, lower_gap / gap
end

_dinf_aspects(dem::AbstractMatrix, cellsize) =
    aspect(dem; method = Horn(), cellsize = abs.(cellsize))
_dinf_aspects(dem, cellsize) = nothing
_dinf_aspect(aspects, dem, cell, cellsize) = aspects[cell]
_dinf_aspect(::Nothing, dem, cell, cellsize) = _localaspect(dem, cell, cellsize)

function _localaspect(dem, cell, cellsize)
    xx = xy = yy = xz = yz = 0.0
    cells = _cells(dem)
    center = dem[cell]
    for neighbor in neighbors(dem, cell)
        neighbor in cells || continue
        distance = celldistance(dem, cell, neighbor; cellsize)
        iszero(distance) && continue
        bearing = _cellbearing(dem, cell, neighbor, cellsize)
        east = distance * sind(bearing)
        north = distance * cosd(bearing)
        elevation = dem[neighbor] - center
        xx += east^2
        xy += east * north
        yy += north^2
        xz += east * elevation
        yz += north * elevation
    end

    determinant = xx * yy - xy^2
    iszero(determinant) && return NaN
    east_gradient = (xz * yy - yz * xy) / determinant
    north_gradient = (yz * xx - xz * xy) / determinant
    iszero(east_gradient) && iszero(north_gradient) && return NaN
    return mod(atand(-east_gradient, -north_gradient), 360)
end

"""
    flowaccumulation(dem::AbstractMatrix, closed=fill!(similar(dem, Bool), false); method=DInf(), cellsize=cellsize(dem))

Computes the flow accumulation of a digital elevation model (DEM) `dem` with an optional `closed` mask and a `method` for flow direction.
Returns the flow accumulation and the flow direction (local drainage direction or ldd).
"""
function flowaccumulation(
    dem,
    closed = fill!(similar(dem, Bool), false);
    method = DInf(),
    cellsize = cellsize(dem),
)
    acc = similar(dem, Float32)
    for cell in _cells(dem)
        acc[cell] = cellarea(dem, cell; cellsize)
    end
    flowaccumulation!(dem, acc, copy(closed); method, cellsize)
end

function flowaccumulation!(
    dem,
    acc::AbstractArray{<:Real},
    closed = fill!(similar(dem, Bool), false);
    method = DInf(),
    cellsize = cellsize(dem),
)
    R = _cells(dem)
    L = similar(dem, Int64)
    L .= LinearIndices(dem)
    first_cell = first(R)
    dir = similar(dem, typeof(first_cell - first_cell))
    fill!(dir, first_cell - first_cell)
    order = ones(Int64, length(closed) - sum(closed))

    open = PriorityQueue{typeof(first_cell), eltype(dem)}()

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        for ncell in neighbors(dem, cell)
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

# Direction grids go back to the caller. The Rasters extensions override this
# to drop the source missingval, which has no representation in the direction
# eltype (Rasters converts it via `typemax` and throws).
_directiongrid(dem, ::Type{C}) where {C <: FlowDirectionConvention} =
    similar(dem, FlowDirection{C, UInt8})

function _accumulate!(::D8, acc, order, dir, R, dem, cellsize)
    for i in reverse(order)
        iszero(dir[i]) && continue
        acc[R[i] + dir[i]] += acc[i]
    end
    output = _directiongrid(dem, LDD)
    output .= FlowDirection{LDD}.(_orient.(dir, Ref(cellsize)))
    return output
end
function _accumulate!(::DInf, acc, order, dir, R, dem, cellsize)
    aspects = _dinf_aspects(dem, cellsize)
    visited = fill!(similar(dem, Bool), false)
    output = _directiongrid(dem, D8D)
    fill!(output, 0)

    for i in reverse(order)
        iszero(dir[i]) && continue
        cell = R[i]
        aspect = _dinf_aspect(aspects, dem, cell, cellsize)

        if !isfinite(aspect)
            acc[cell + dir[i]] += acc[i]
            output[i] = _directioncode(dem, dir[i], cellsize)
            visited[i] = true
            continue
        end

        flowneighbors = _flowneighbors(dem, cell, aspect, cellsize)
        if isnothing(flowneighbors)
            acc[cell + dir[i]] += acc[i]
            output[i] = _directioncode(dem, dir[i], cellsize)
            visited[i] = true
            continue
        end
        acell, bcell, aa, bb = flowneighbors
        a, b = acell - cell, bcell - cell

        # Depression
        if (a != dir[i] && b != dir[i])
            acc[cell + dir[i]] += acc[i]
            output[i] = _directioncode(dem, dir[i], cellsize)
            visited[i] = true
            continue
        end

        # Scale flows correctly at the edges
        if !(acell in R) || visited[acell]
            aa = 0
            bb = 1
        end
        if !(bcell in R) || visited[bcell]
            aa = 1
            bb = 0
        end

        dirs = zero(UInt8)
        if acell in R && aa > 0
            acc[acell] += acc[i] * aa
            dirs |= _directioncode(dem, a, cellsize)
        end
        if bcell in R && bb > 0
            acc[bcell] += acc[i] * bb
            dirs |= _directioncode(dem, b, cellsize)
        end
        output[i] = dirs
        if acell in R && bcell in R && visited[acell] && visited[bcell]
            error()
        end
        visited[i] = true
    end
    return output
end

function _accumulate!(fd8::FD8, acc, order, dir, R, dem::AbstractMatrix, cellsize)
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
    output = _directiongrid(dem, D8D)
    fill!(output, 0)

    weights = zeros(size(contour_lengths))
    Σw = 0.0
    for i in reverse(order)
        dir[i] == CartesianIndex(0, 0) && continue
        fill!(weights, 0)
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
    end
    return output
end

function _accumulate!(fd8::FD8, acc, order, dir, R, dem, cellsize)
    visited = fill!(similar(dem, Bool), false)
    output = _directiongrid(dem, D8D)
    fill!(output, 0)

    for i in reverse(order)
        iszero(dir[i]) && continue
        cell = R[i]
        total_weight = 0.0
        downstream = cell + dir[i]
        downstream_weight = 0.0
        for neighbor in neighbors(dem, cell)
            neighbor in R || continue
            weight = _fd8weight(fd8, dem, cell, neighbor, R, visited, cellsize)
            total_weight += weight
            neighbor == downstream && (downstream_weight = weight)
        end

        if iszero(total_weight) || iszero(downstream_weight)
            acc[downstream] += acc[i]
            output[i] = _directioncode(dem, dir[i], cellsize)
            visited[i] = true
            continue
        end

        directions = zero(UInt8)
        for neighbor in neighbors(dem, cell)
            neighbor in R || continue
            weight = _fd8weight(fd8, dem, cell, neighbor, R, visited, cellsize)
            iszero(weight) && continue
            acc[neighbor] += acc[i] * (weight / total_weight)
            directions |= _directioncode(dem, neighbor - cell, cellsize)
        end
        output[i] = directions
        visited[i] = true
    end
    return output
end

function _fd8weight(fd8, dem, cell, neighbor, cells, visited, cellsize)
    visited[neighbor] && return 0.0
    difference = dem[cell] - dem[neighbor]
    (difference < 0 || isnan(difference)) && return 0.0

    distance = celldistance(dem, cell, neighbor; cellsize)
    iszero(distance) && return 0.0
    bearing = _cellbearing(dem, cell, neighbor, cellsize)
    previous_gap = next_gap = Inf
    for other in neighbors(dem, cell)
        (other == neighbor || !(other in cells)) && continue
        other_bearing = _cellbearing(dem, cell, other, cellsize)
        candidate_previous_gap = mod(bearing - other_bearing, 360)
        candidate_next_gap = mod(other_bearing - bearing, 360)
        iszero(candidate_previous_gap) ||
            (previous_gap = min(previous_gap, candidate_previous_gap))
        iszero(candidate_next_gap) || (next_gap = min(next_gap, candidate_next_gap))
    end
    (!isfinite(previous_gap) || !isfinite(next_gap)) && return 0.0

    contour_length =
        distance / 2 * (tand(previous_gap / 2) + tand(next_gap / 2))
    contour_length > 0 || return 0.0
    return (difference / distance * contour_length)^fd8.p
end

"""
    topographic_wetness_index(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))
    topographic_wetness_index(dem; method=D8(), cellsize=cellsize(dem))

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
function topographic_wetness_index(
    dem;
    method = D8(),
    cellsize = cellsize(dem),
)
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc / tand(s))
end
@deprecate TWI topographic_wetness_index

"""
    stream_power_index(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))
    stream_power_index(dem; method=D8(), cellsize=cellsize(dem))

Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function stream_power_index(dem::AbstractMatrix; method = DInf(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end
function stream_power_index(dem; method = D8(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end
@deprecate SPI stream_power_index

"""
    height_above_nearest_drainage(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem), threshold=100)

Compute Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem` 
with an optional `method` for flow direction, a `cellsize`, and a flow accumulation `threshold` for stream definition.
"""
function height_above_nearest_drainage(
    dem;
    method = D8(),
    cellsize = cellsize(dem),
    threshold = 100,
)
    R = _cells(dem)
    L = similar(dem, Int64)
    L .= LinearIndices(dem)
    first_cell = first(R)
    dir = similar(dem, typeof(first_cell - first_cell))
    fill!(dir, first_cell - first_cell)
    closed = similar(dem, Bool)
    closed .= false
    order = ones(Int64, length(closed) - sum(closed))

    output = zero(dem)
    acc = similar(dem, Float32)
    for cell in R
        acc[cell] = cellarea(dem, cell; cellsize)
    end

    open = PriorityQueue{typeof(first_cell), eltype(dem)}()

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        for ncell in neighbors(dem, cell)
            ncell in R || continue
            # skip visited and center cells
            (closed[ncell] || ncell == cell) && continue

            closed[ncell] = true
            dir[ncell] = cell - ncell

            enqueue!(open, ncell, dem[ncell])
        end
    end

    flowdirs = _accumulate!(method, acc, order, dir, R, dem, cellsize)
    stream_mask = acc .>= threshold
    _burn_streams!(stream_mask, order, dir, R)
    _fill_flow_depressions!(acc, order, flowdirs, R, dem, cellsize)
    _hand!(output, order, flowdirs, R, acc, stream_mask, cellsize)
    return output
end

"""
    height_above_nearest_drainage(dem::AbstractMatrix, stream_mask::AbstractMatrix{Bool})

Computes the Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem` 
given a stream definition as a boolean `stream_mask`.
"""
function height_above_nearest_drainage(
    dem,
    stream_mask::AbstractArray{Bool},
)
    R = _cells(dem)
    L = LinearIndices(dem)
    first_cell = first(R)
    dir = fill(first_cell - first_cell, size(dem))
    closed = fill!(similar(dem, Bool), false)
    order = ones(Int64, length(closed) - sum(closed))

    output = zero(dem)

    open = PriorityQueue{typeof(first_cell), eltype(dem)}()

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        for ncell in neighbors(dem, cell)
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
            # to the downstream cell's HAND value.
            output[i] = output[R[i] + dir[i]] + dem[i] - dem[R[i] + dir[i]]
        else
            output[i] = NaN
        end
    end
    for i in eachindex(output)
        output[i] = max(zero(output[i]), output[i])
    end
    return output
end

function _burn_streams!(stream_mask, order, dir, R)
    for i in reverse(order)
        stream_mask[i] || continue
        iszero(dir[i]) && continue
        stream_mask[R[i] + dir[i]] = true
    end
    return stream_mask
end

function _fill_flow_depressions!(filled, order, flowdirs, R, dem, cellsize)
    RelativeIndex = typeof(first(R) - first(R))
    fill!(filled, NaN)
    for i in order
        cell = R[i]
        elevation = dem[i]
        if !isfinite(elevation)
            filled[i] = NaN
            continue
        end
        if !ispit(flowdirs[i])
            for relative in _relativeindices(RelativeIndex, flowdirs[i], cell, cellsize)
                downstream = cell + relative
                downstream in R || continue
                isfinite(filled[downstream]) || continue
                elevation = max(elevation, filled[downstream])
            end
        end
        filled[i] = elevation
    end
    return filled
end

function _hand!(output, order, flowdirs, R, filled, stream_mask, cellsize)
    RelativeIndex = typeof(first(R) - first(R))
    for i in order
        cell = R[i]
        if stream_mask[i] || ispit(flowdirs[i])
            output[i] = 0.0
        elseif !isfinite(filled[i])
            output[i] = NaN
        else
            found = false
            best = zero(eltype(output))
            for relative in _relativeindices(RelativeIndex, flowdirs[i], cell, cellsize)
                downstream = cell + relative
                downstream in R || continue
                if !isfinite(filled[downstream]) || !isfinite(output[downstream])
                    continue
                end
                candidate =
                    output[downstream] + filled[i] - filled[downstream]
                if !found || candidate < best
                    best = candidate
                    found = true
                end
            end
            output[i] = found ? best : NaN
        end
    end
    for i in eachindex(output)
        output[i] = max(zero(output[i]), output[i])
    end
    return output
end

"""
    depression_depth(dem; filled=filldepressions(dem))

Computes the depth of each cell below the filled surface.

Returns the difference between the depression-filled DEM and the original DEM,
representing how deep each cell sits within a depression. Cells not in
depressions will have a depth of zero.

This is useful for identifying potential cold air pooling zones and water
retention areas.
"""
function depression_depth(dem; filled = filldepressions(dem))
    filled .- dem
end

"""
    depression_volume(dem::AbstractMatrix; filled=filldepressions(dem), cellsize=cellsize(dem))
    depression_volume(dem; filled=filldepressions(dem), cellsize=cellsize(dem))

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
function depression_volume(
    dem;
    filled = filldepressions(dem),
    cellsize = cellsize(dem),
)
    depth = depression_depth(dem; filled)
    volume = 0.0
    for cell in _cells(dem)
        isfinite(depth[cell]) || continue
        volume += depth[cell] * cellarea(dem, cell; cellsize)
    end
    return volume
end

"""
    drainage_potential(dem::AbstractMatrix; method=DInf(), cellsize=cellsize(dem))
    drainage_potential(dem; method=D8(), cellsize=cellsize(dem))

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
function drainage_potential(dem; method = D8(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. sind(s) / log1p(acc)
end
