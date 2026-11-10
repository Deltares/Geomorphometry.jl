
function outlets(A::AbstractMatrix)
    CI = CartesianIndices(A)
    indices = Vector{CartesianIndex}()
    append!(indices, first(eachrow(CI)))
    append!(indices, last(eachrow(CI)))
    append!(indices, first(eachcol(CI)))
    append!(indices, last(eachcol(CI)))
    unique(indices)
end

cellarea(dem::AbstractMatrix, cell; cellsize = cellsize(dem)) =
    abs(cellsize[1] * cellsize[2])

abstract type DepressionFillingMethod end

"""
    PriorityFlood()

Improved Priority-Flood, Algorithm 2 of
[Barnes et al. (2014)](@cite barnesPriorityFloodOptimalDepressionFilling2014). Cells are
flooded inwards from the outlets through a priority queue, and cells found below the
current flood level are filled through a plain queue. Runs in `O(m log m)` time, where
`m` is the number of cells that reach the priority queue.

The flood itself is that of [Wang and Liu (2006)](@cite wangEfficientMethodIdentifying2006),
who derived it as a least-cost search whose cost is a cell's *spill elevation*, the lowest
level it must be raised to for water to reach an outlet. The plain queue is what Barnes et
al. add, keeping cells below the flood level out of the priority queue, for run times 16%
shorter.
"""
struct PriorityFlood <: DepressionFillingMethod end

"""
    Zhou()

Priority-Flood variant of [Zhou et al. (2016)](@cite zhouEfficientVariantPriorityFlood2016).
Terrain that rises away from an already flooded cell needs no priority ordering, so slope
cells are region-grown on a plain queue and only *potential spill cells*, where the terrain
turns back downhill, reach the priority queue. Runs in `O(N + m log m)` time with `m` much
smaller than the number of cells `N`, and reported 44.6% faster than
[`PriorityFlood`](@ref). Measured here it is worth 1.1-1.4x, most on continuous elevations
with few ties and least on tiles with much nodata.
"""
struct Zhou <: DepressionFillingMethod end

"""
    filldepressions(dem, mask=fill!(similar(dem, Bool), false); method=PriorityFlood(), flats=true)

Fill the depressions of a digital elevation model, so that every cell drains to an outlet.

# Arguments
- `dem`: An indexed digital elevation model.
- `mask`: An optional boolean array representing the mask. Cells with `true` values are treated as already filled (queued) and excluded from the computation. Defaults to all `false`.
- `method::DepressionFillingMethod=PriorityFlood()`: Algorithm used, either
  [`PriorityFlood`](@ref) or [`Zhou`](@ref). The two produce the same filled surface for a
  DEM of finite elevations. They do **not** agree once the DEM carries nodata: neither
  algorithm defines an ordering for `NaN`, so which cells a flood reaches, and in what
  order, is left to the comparison rules rather than decided. Mask nodata out beforehand if
  the choice of method must not change the result.
- `flats`: Also give the resulting flats a drainage gradient with [`resolveflats`](@ref),
  so that every cell of the output has a strictly lower neighbour and a defined flow
  direction. Pass `false` for the filled surface on its own. Defaults to `true` for a
  floating point `dem` and `false` otherwise, as integer elevations cannot express a
  gradient of single float steps. Cells excluded by `mask` keep their elevation.
"""
function filldepressions(
    dem,
    mask = fill!(similar(dem, Bool), false);
    method::DepressionFillingMethod = PriorityFlood(),
    flats::Bool = _canresolveflats(dem),
)
    filldepressions!(copy(dem), mask; method, flats)
end

# Resolving flats nudges elevations by single float steps, which an integer DEM cannot
# represent, so it is only on by default where it can work.
_canresolveflats(dem) = eltype(dem) <: AbstractFloat

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

function filldepressions!(
    dem,
    queued = fill!(similar(dem, Bool), false);
    method::DepressionFillingMethod = PriorityFlood(),
    flats::Bool = _canresolveflats(dem),
)
    # The flood writes into `queued`, so the caller's exclusions are kept to hand back to
    # `resolveflats!`, which must leave those cells alone too.
    excluded = flats ? copy(queued) : nothing
    _filldepressions!(method, dem, queued)
    flats && resolveflats!(dem, excluded)
    return dem
end

function _filldepressions!(::PriorityFlood, dem, queued)
    R = eachcell(dem)
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
            queued[ncell] && continue
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

# Zhou et al. (2016). Three queues: `open` holds only potential spill cells, while `trace`
# and `depression` region-grow slope and depression cells on plain queues. A cell reached
# by walking uphill from a processed cell drains back the way it came, so its elevation
# does not depend on visit order and it never needs the priority queue.
function _filldepressions!(::Zhou, dem, processed)
    first_cell = first(eachcell(dem))
    Cell = typeof(first_cell)
    open = PriorityQueue{Cell, eltype(dem)}()
    trace = DataStructures.Queue{Cell}()
    depression = DataStructures.Queue{Cell}()

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        processed[cell] = true
    end

    @inbounds while !isempty(open)
        cell = dequeue!(open)
        spill = dem[cell]
        for ncell in neighbors(dem, cell)
            processed[ncell] && continue
            processed[ncell] = true
            if dem[ncell] <= spill
                dem[ncell] = spill
                DataStructures.enqueue!(depression, ncell)
                _processdepression!(dem, processed, depression, trace)
            else
                DataStructures.enqueue!(trace, ncell)
            end
            _processtrace!(dem, processed, trace, open)
        end
    end
    return dem
end

"""Region-grow a depression at the spill elevation, handing higher cells to `trace`."""
function _processdepression!(dem, processed, depression, trace)
    @inbounds while !isempty(depression)
        cell = DataStructures.dequeue!(depression)
        spill = dem[cell]
        for ncell in neighbors(dem, cell)
            processed[ncell] && continue
            processed[ncell] = true
            if dem[ncell] > spill
                DataStructures.enqueue!(trace, ncell)
            else
                dem[ncell] = spill
                DataStructures.enqueue!(depression, ncell)
            end
        end
    end
    return dem
end

"""
Region-grow slope cells uphill on a plain queue. A cell is pushed onto the priority
queue only when it borders unprocessed terrain at or below its own elevation and no
processed neighbour of that terrain is lower, which makes it a potential spill cell.
"""
function _processtrace!(dem, processed, trace, open)
    @inbounds while !isempty(trace)
        cell = DataStructures.dequeue!(trace)
        spill = dem[cell]
        queued = false
        for ncell in neighbors(dem, cell)
            processed[ncell] && continue
            nspill = dem[ncell]
            if nspill <= spill
                if !queued && _isspillcell(dem, processed, ncell, nspill)
                    enqueue!(open, cell, spill)
                    queued = true
                end
                continue
            end
            processed[ncell] = true
            DataStructures.enqueue!(trace, ncell)
        end
    end
    return dem
end

@inline function _isspillcell(dem, processed, ncell, nspill)
    @inbounds for jcell in neighbors(dem, ncell)
        processed[jcell] && dem[jcell] < nspill && return false
    end
    return true
end

# Non-diagonal neighbors come first: Algorithm 4 of
# [barnesPriorityFloodOptimalDepressionFilling2014](@cite) processes them first on its
# line 13, as they have the greatest centre-to-centre slopes and so take priority when
# several neighbors sit at the same elevation.
const nbs =
    CartesianIndex.([(-1, 0), (0, -1), (0, 1), (1, 0), (-1, -1), (-1, 1), (1, -1), (1, 1)])

"""
    neighbors(dem, cell)

Iterate over the neighbors of `cell` that fall within `dem`. Matrices use the
Moore neighborhood. Implementations must yield only in-domain indices: kernels
index `dem` with the result directly, some under `@inbounds`.
"""
neighbors(dem, cell::CartesianIndex{2}) =
    Iterators.filter(in(CartesianIndices(dem)), cell + nb for nb in nbs)

function watersheds(dem, queued = fill!(similar(dem, Bool), false))
    first_cell = first(eachcell(dem))
    open = PriorityQueue{typeof(first_cell), eltype(dem)}()
    pit = DataStructures.Queue{typeof(first_cell)}()
    labels = fill!(similar(dem, Int), 0)
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
            queued[ncell] && continue
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

# The eight triangular facets of [Tarboton (1997)](@cite tarbotonNewMethodDetermination1997),
# each pairing an orthogonal neighbour with one of the two diagonals beside it.
const _facets = (
    (CartesianIndex(1, 0), CartesianIndex(1, 1)),
    (CartesianIndex(0, 1), CartesianIndex(1, 1)),
    (CartesianIndex(0, 1), CartesianIndex(-1, 1)),
    (CartesianIndex(-1, 0), CartesianIndex(-1, 1)),
    (CartesianIndex(-1, 0), CartesianIndex(-1, -1)),
    (CartesianIndex(0, -1), CartesianIndex(-1, -1)),
    (CartesianIndex(0, -1), CartesianIndex(1, -1)),
    (CartesianIndex(1, 0), CartesianIndex(1, -1)),
)

"""
    _dinfreceivers(dem, cell, state, cellsize)

The two neighbours [`DInf`](@ref) splits `cell`'s flow between, and their fractions, or
`nothing` when it has nowhere to send it.

On a raster this is the facet method of
[Tarboton (1997)](@cite tarbotonNewMethodDetermination1997): the steepest of the eight
triangular facets around the cell fixes a flow angle, which is shared between the two
neighbours bounding that facet. Other grids have no facet geometry, so they fall back to
bracketing the local aspect.
"""
function _dinfreceivers(dem::AbstractMatrix, cell, state, cellsize)
    R = CartesianIndices(dem)
    center = dem[cell]
    isfinite(center) || return nothing

    steepest = 0.0
    near = far = cell
    slope1 = slope2 = length1 = length2 = 0.0
    onedge = 0        # -1 towards the orthogonal neighbour, +1 towards the diagonal

    @inbounds for (offset1, offset2) in _facets
        first_cell, second_cell = cell + offset1, cell + offset2
        (first_cell in R && second_cell in R) || continue
        z1, z2 = dem[first_cell], dem[second_cell]
        (isfinite(z1) && isfinite(z2)) || continue
        d1 = celldistance(dem, cell, first_cell; cellsize)
        d2 = celldistance(dem, first_cell, second_cell; cellsize)
        (d1 > 0 && d2 > 0) || continue

        s1 = (center - z1) / d1
        s2 = (z1 - z2) / d2
        # Where the steepest direction of the facet falls outside it, the maximum sits on
        # one of its edges. Both tests are comparisons: no trigonometry is needed to find
        # out which case holds, only to share the flow once the winner is known.
        if s2 <= 0
            slope, side = s1, -1
        elseif s2 * d1 >= s1 * d2
            slope, side = (center - z2) / hypot(d1, d2), 1
        else
            slope, side = hypot(s1, s2), 0
        end

        if slope > steepest
            steepest = slope
            near, far = first_cell, second_cell
            slope1, slope2, length1, length2 = s1, s2, d1, d2
            onedge = side
        end
    end

    steepest > 0 || return nothing
    fraction = onedge < 0 ? 0.0 :
        onedge > 0 ? 1.0 : atan(slope2, slope1) / atan(length2, length1)
    return near, far, 1 - fraction, fraction
end

function _dinfreceivers(dem, cell, state, cellsize)
    aspect = _dinf_aspect(state, dem, cell, cellsize)
    # Undefined over a flat or beside nodata, and bracketing a NaN would pick the cell
    # itself as both of its receivers.
    isfinite(aspect) || return nothing
    return _flowneighbors(dem, cell, aspect, cellsize)
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

# Rasters take the facet route, which needs no precomputation; other grids derive an aspect
# per cell. Either way there is nothing to prepare for the whole grid.
_dinf_aspects(dem, cellsize) = nothing
_dinf_aspect(aspects, dem, cell, cellsize) = aspects[cell]
_dinf_aspect(::Nothing, dem, cell, cellsize) = _localaspect(dem, cell, cellsize)

function _localaspect(dem, cell, cellsize)
    east_gradient, north_gradient = _symmetric_gradient(
        SymmetricGradient(),
        dem,
        cell,
        dem[cell],
        _neighborvalues(dem, cell),
        cellsize,
    )
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
    for cell in eachcell(dem)
        acc[cell] = cellarea(dem, cell; cellsize)
    end
    flowaccumulation!(dem, acc, copy(closed); method, cellsize)
end

"""
    _carve(dem, closed, cellsize)

Flood `dem` inwards from its outlets and return the order cells were reached in, together
with the direction each one drains to. `dem` is not modified: flow is routed through
depressions rather than over them.

Algorithm 4 of
[Barnes et al. (2014)](@cite barnesPriorityFloodOptimalDepressionFilling2014), after the
least-cost path search of [Metz et al. (2011)](@cite metzEfficientExtractionDrainage2011).
A cell settles its direction when it is popped, taking the steepest of the neighbours
already popped: the flood goes in elevation order, so those are exactly the ones below it.

Cells already marked in `closed` are left out of the order, and `closed` is written to as
the flood proceeds.
"""
function _carve(dem, closed, cellsize)
    R = eachcell(dem)
    L = similar(dem, Int64)
    L .= LinearIndices(dem)
    first_cell = first(R)
    dir = _alloc(dem, typeof(first_cell - first_cell))
    fill!(dir, first_cell - first_cell)
    order = ones(Int64, length(closed) - sum(closed))

    open = PriorityQueue{typeof(first_cell), eltype(dem)}()
    drained = fill!(similar(dem, Bool), false)  # popped, and so downstream of what follows
    rim = fill!(similar(dem, Bool), false)      # flow leaves the grid here, not a cell

    @inbounds for cell in outlets(dem)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
        rim[cell] = true
    end
    i = 1
    @inbounds while !isempty(open)
        cell = dequeue!(open)
        order[i] = L[cell]
        i += 1
        drained[cell] = true

        z = dem[cell]
        # Nothing seen yet, so the first neighbour already popped wins outright and a
        # level one stands until something lower turns up, which is how a cell on flat
        # ground still drains.
        steepest = typemin(float(eltype(dem)))
        for ncell in neighbors(dem, cell)
            if drained[ncell]
                # Already popped, so it is downstream: a candidate receiver for this cell.
                rim[cell] && continue
                nz = dem[ncell]
                # A -Inf neighbour is the lowest thing there is and a legitimate receiver,
                # which the manual relies on to drain into sea. NaN and +Inf are not.
                isnan(nz) && continue
                (isinf(nz) && !signbit(nz)) && continue
                distance = celldistance(dem, cell, ncell; cellsize)
                iszero(distance) && continue
                gradient = (z - nz) / distance
                if gradient > steepest
                    steepest = gradient
                    dir[cell] = ncell - cell
                end
            elseif !closed[ncell]
                closed[ncell] = true
                enqueue!(open, ncell, dem[ncell])
            end
        end
    end
    return order, dir
end

function flowaccumulation!(
    dem,
    acc::AbstractArray{<:Real},
    closed = fill!(similar(dem, Bool), false);
    method = DInf(),
    cellsize = cellsize(dem),
)
    R = eachcell(dem)
    order, dir = _carve(dem, closed, cellsize)

    dirs = _accumulate!(method, acc, order, dir, R, dem, cellsize)
    return acc, dirs
end

"""
Allocate a direction grid shaped like `dem`. A hook, because a direction eltype has no
representation for a source `missingval` and some grid types must drop theirs.
"""
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

        flowneighbors = _dinfreceivers(dem, cell, aspects, cellsize)
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

        # Do not route flow back into cells that have already drained
        if visited[acell]
            aa = 0
            bb = 1
        end
        if visited[bcell]
            aa = 1
            bb = 0
        end

        dirs = zero(UInt8)
        if aa > 0
            acc[acell] += acc[i] * aa
            dirs |= _directioncode(dem, a, cellsize)
        end
        if bb > 0
            acc[bcell] += acc[i] * bb
            dirs |= _directioncode(dem, b, cellsize)
        end
        output[i] = dirs
        if visited[acell] && visited[bcell]
            error()
        end
        visited[i] = true
    end
    return output
end

"""
Rasters keep the fixed distance and contour tables of
[Quinn et al. (1991)](@cite quinnPredictionHillslopeFlow1991), which depend only on
`cellsize`. Grids with no such fixed neighbourhood derive both per cell instead, in
[`_fd8weights!`](@ref).
"""
function _accumulate!(fd8::FD8, acc, order, dir, R::CartesianIndices{2}, dem, cellsize)
    # Contour lengths come from `_rastercontour`, which `flowdirection` uses too.
    offsets = CartesianIndices((-1:1, -1:1))
    contour_lengths = map(o -> iszero(o) ? 0.0 : _rastercontour(o, cellsize), offsets)

    δx, δy = abs.(cellsize)
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
    scratch = _fd8scratch(dem)

    for i in reverse(order)
        iszero(dir[i]) && continue
        cell = R[i]
        downstream = cell + dir[i]
        total_weight = _fd8weights!(scratch, fd8, dem, cell, visited, cellsize)
        downstream_weight = 0.0
        for (index, neighbor) in enumerate(scratch.neighbors)
            neighbor == downstream && (downstream_weight = scratch.weights[index])
        end

        if iszero(total_weight) || iszero(downstream_weight)
            acc[downstream] += acc[i]
            output[i] = _directioncode(dem, dir[i], cellsize)
            visited[i] = true
            continue
        end

        directions = zero(UInt8)
        for (index, neighbor) in enumerate(scratch.neighbors)
            weight = scratch.weights[index]
            iszero(weight) && continue
            acc[neighbor] += acc[i] * (weight / total_weight)
            directions |= _directioncode(dem, neighbor - cell, cellsize)
        end
        output[i] = directions
        visited[i] = true
    end
    return output
end

"""
Scratch space for [`_fd8weights!`](@ref), reused across cells.

Bearings are computed once per cell rather than once per neighbour pair, and the contour
lengths they feed depend only on the *shape* of the neighbourhood, not on elevation, so
they are recomputed only when that shape changes — on a raster, only at the edges.
"""
mutable struct FD8Scratch{C, O, T}
    neighbors::Vector{C}
    offsets::Vector{O}
    bearings::Vector{T}
    weights::Vector{T}
    contours::Vector{T}
    distances::Vector{T}
    shape::Vector{O}
    bitweights::Vector{T}   # a cell's shares, gathered against their direction bit
end

"""
    _rastercontour(offset, cellsize)

Contour length a step of `offset` crosses on a square grid, after
[Quinn et al. (1991)](@cite quinnPredictionHillslopeFlow1991).

The eight Moore neighbours are not a Voronoi neighbourhood — a diagonal shares only a
corner, so there is no face there to measure — and the angular construction in
[`_contourlengths!`](@ref) degenerates on them, giving every direction the same weight. This
is the fix-up for that, and both [`flowdirection`](@ref) and [`flowaccumulation`](@ref) take
their raster contour lengths from here so the two cannot disagree.
"""
function _rastercontour(offset, cellsize)
    δx, δy = abs.(cellsize)
    (!iszero(offset[1]) && !iszero(offset[2])) && return hypot(δx, δy) / 4
    # Half the face the step crosses, which is perpendicular to it: a step along axis 1
    # travels δx but leaves through a face δy wide. `dists` beside the caller asks the
    # other question, how far the step goes, and so pairs each axis with its own spacing.
    return iszero(offset[2]) ? δy / 2 : δx / 2
end

function _fd8scratch(dem)
    cell = first(eachcell(dem))
    C, O = typeof(cell), typeof(cell - cell)
    FD8Scratch(C[], O[], Float64[], Float64[], Float64[], Float64[], O[], zeros(Float64, 8))
end

"""
Recompute the contour length and distance of each neighbour for the neighbourhood shape now
in `scratch.offsets`. Each neighbour is given the Voronoi edge it would share with `cell`,
derived from the angular gaps to its neighbours on either side, which is exact wherever the
neighbourhood is a Voronoi one. Rasters are the exception and have their own method.
"""
function _contourlengths!(scratch::FD8Scratch, dem, cell, cellsize)
    n = length(scratch.offsets)
    resize!(scratch.contours, n)
    resize!(scratch.distances, n)
    bearings = scratch.bearings
    resize!(bearings, n)
    @inbounds for (index, neighbor) in enumerate(scratch.neighbors)
        bearings[index] = _cellbearing(dem, cell, neighbor, cellsize)
        scratch.distances[index] = celldistance(dem, cell, neighbor; cellsize)
    end

    @inbounds for index in 1:n
        bearing = bearings[index]
        previous_gap = next_gap = Inf
        for other in 1:n
            other == index && continue
            candidate_previous_gap = mod(bearing - bearings[other], 360)
            candidate_next_gap = mod(bearings[other] - bearing, 360)
            iszero(candidate_previous_gap) ||
                (previous_gap = min(previous_gap, candidate_previous_gap))
            iszero(candidate_next_gap) || (next_gap = min(next_gap, candidate_next_gap))
        end
        # A gap wider than a right angle means neighbours are missing, at the edge of the
        # DEM or beside nodata. Clamp it: the formula otherwise runs the contour out to
        # tand(90) = Inf, and an infinite weight takes the whole cell's flow.
        #
        # This models the boundary as a flat face perpendicular to the centre-to-centre
        # line at half the distance, which is a Voronoi edge. So it is exact wherever
        # the neighbourhood is a Voronoi one: d for the four faces of a 4-connected
        # square, d/√3 for the six edges of a hexagon, every contour equal on a
        # hexagonal DGGS.
        #
        # TODO The 8-connected Moore neighbourhood is the one case that is not a Voronoi
        # neighbourhood — diagonal neighbours share only a corner, no edge — so there is
        # no face here to measure and the formula returns a degenerate answer: equal
        # gaps force L = d * tand(gap / 2), exactly proportional to d, and the 1/√2 a
        # diagonal loses on slope it wins straight back on contour. Every direction then
        # weighs 0.414Δz, ranking neighbours by drop rather than gradient, with `p`
        # concentrating on the largest drop rather than the steepest slope.
        #
        # The missing term is the projection: flux through a boundary element is
        # S * w * cos(angle between the flow and that element's normal), and summing raw
        # arc length assumes the boundary faces the flow squarely. It does on a Voronoi
        # edge, which is why this is exact there and why the term is invisible until a
        # corner shows up. Projecting turns the diagonal's 0.586 into 0.414 and restores
        # a √2 ratio; Quinn's raster table splits each face 50/25/25 instead and
        # projects to 0.5 and 0.354, a ratio of 2. Hence the separate raster method.
        scratch.contours[index] =
            (isfinite(previous_gap) && isfinite(next_gap)) ?
            scratch.distances[index] / 2 *
            (tand(min(previous_gap, 90) / 2) + tand(min(next_gap, 90) / 2)) : 0.0
    end

    resize!(scratch.shape, n)
    copyto!(scratch.shape, scratch.offsets)
    return scratch
end

"""Rasters take their contour lengths from [`_rastercontour`](@ref)."""
function _contourlengths!(
    scratch::FD8Scratch{<:CartesianIndex{2}}, dem, cell, cellsize
)
    n = length(scratch.offsets)
    resize!(scratch.contours, n)
    resize!(scratch.distances, n)
    @inbounds for (index, offset) in enumerate(scratch.offsets)
        scratch.distances[index] =
            celldistance(dem, cell, scratch.neighbors[index]; cellsize)
        scratch.contours[index] = _rastercontour(offset, cellsize)
    end
    resize!(scratch.shape, n)
    copyto!(scratch.shape, scratch.offsets)
    return scratch
end

"""
    _fd8weights!(scratch, fd8, dem, cell, visited, cellsize)

Fill `scratch` with the neighbours of `cell` and their unnormalised FD8 weights, and return
the total. A neighbour weighs the downslope gradient towards it times the length of the
cell's contour facing it, raised to `fd8.p`, after
[Quinn et al. (1991)](@cite quinnPredictionHillslopeFlow1991). Neighbours that are not
lower, or that `visited` marks as already processed, weigh nothing.
"""
function _fd8weights!(scratch::FD8Scratch, fd8, dem, cell, visited, cellsize)
    empty!(scratch.neighbors)
    empty!(scratch.offsets)
    @inbounds for neighbor in neighbors(dem, cell)
        push!(scratch.neighbors, neighbor)
        push!(scratch.offsets, neighbor - cell)
    end
    scratch.offsets == scratch.shape || _contourlengths!(scratch, dem, cell, cellsize)
    resize!(scratch.weights, length(scratch.neighbors))

    center = dem[cell]
    total = 0.0
    @inbounds for (index, neighbor) in enumerate(scratch.neighbors)
        weight = 0.0
        difference = center - dem[neighbor]
        # No drop, no flow. Skipping here also keeps a zero difference away from a contour
        # length of zero, which a truncated neighbourhood can produce.
        if !visited[neighbor] && difference > 0 && !isnan(difference)
            distance = scratch.distances[index]
            contour_length = scratch.contours[index]
            if distance > 0 && contour_length > 0
                weight = (difference / distance * contour_length)^fd8.p
                total += weight
            end
        end
        scratch.weights[index] = weight
    end
    return total
end

"""
    topographic_wetness_index(dem; method=DInf(), cellsize=cellsize(dem))

Computes the Topographic Wetness Index (TWI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function topographic_wetness_index(
    dem;
    method = DInf(),
    cellsize = cellsize(dem),
)
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc / tand(s))
end
@deprecate TWI topographic_wetness_index

"""
    stream_power_index(dem; method=DInf(), cellsize=cellsize(dem))

Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function stream_power_index(dem; method = DInf(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end
@deprecate SPI stream_power_index

"""
    height_above_nearest_drainage(dem; method=DInf(), cellsize=cellsize(dem), threshold=100)

Compute Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem`
with an optional `method` for flow direction, a `cellsize`, and a flow accumulation `threshold` for stream definition.
"""
function height_above_nearest_drainage(
    dem;
    method = DInf(),
    cellsize = cellsize(dem),
    threshold = 100,
)
    R = eachcell(dem)
    closed = fill!(similar(dem, Bool), false)
    output = zero(dem)
    acc = similar(dem, Float32)
    for cell in R
        acc[cell] = cellarea(dem, cell; cellsize)
    end

    order, dir = _carve(dem, closed, cellsize)

    flowdirs = _accumulate!(method, acc, order, dir, R, dem, cellsize)
    stream_mask = acc .>= threshold
    _burn_streams!(stream_mask, order, dir, R)
    # At least the precision of the DEM, and floating point to hold NaN
    filled = similar(dem, promote_type(Float32, eltype(dem)))
    _fill_flow_depressions!(filled, order, flowdirs, R, dem, cellsize)
    _hand!(output, order, flowdirs, R, filled, stream_mask, cellsize)
    return output
end

"""
    height_above_nearest_drainage(dem, stream_mask::AbstractArray{Bool})

Computes the Height Above Nearest Drainage (HAND, [nobreHeightNearestDrainage2011](@cite)) of a digital elevation model (DEM) `dem`
given a stream definition as a boolean `stream_mask`.
"""
function height_above_nearest_drainage(
    dem,
    stream_mask::AbstractArray{Bool};
    cellsize = cellsize(dem),
)
    R = eachcell(dem)
    closed = fill!(similar(dem, Bool), false)
    output = zero(dem)

    order, dir = _carve(dem, closed, cellsize)

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
    depression_depth(dem; filled=filldepressions(dem; flats=false))

Computes the depth of each cell below the filled surface.

Returns the difference between the depression-filled DEM and the original DEM,
representing how deep each cell sits within a depression. Cells not in
depressions will have a depth of zero.

This is useful for identifying potential cold air pooling zones and water
retention areas.
"""
function depression_depth(dem; filled = filldepressions(dem; flats=false))
    filled .- dem
end

"""
    depression_volume(dem::AbstractMatrix; filled=filldepressions(dem; flats=false), cellsize=cellsize(dem))
    depression_volume(dem; filled=filldepressions(dem; flats=false), cellsize=cellsize(dem))

Computes the total volume of all depressions/basins in the DEM.

Returns the sum of depression depths multiplied by cell area.
"""
function depression_volume(
    dem::AbstractMatrix;
    filled = filldepressions(dem; flats = false),
    cellsize = cellsize(dem),
)
    sum(filter(isfinite, depression_depth(dem; filled))) * prod(abs.(cellsize))
end
function depression_volume(
    dem;
    filled = filldepressions(dem; flats = false),
    cellsize = cellsize(dem),
)
    depth = depression_depth(dem; filled)
    volume = 0.0
    for cell in eachcell(dem)
        isfinite(depth[cell]) || continue
        volume += depth[cell] * cellarea(dem, cell; cellsize)
    end
    return volume
end

"""
    drainage_potential(dem; method=DInf(), cellsize=cellsize(dem))

Computes a drainage potential index indicating how well each cell drains.

High values indicate good drainage (steep slopes with low upstream accumulation).
Low values indicate poor drainage (flat areas that accumulate flow from upslope).

This is useful as a proxy for cold air drainage - cells with high drainage potential
will shed cold air downslope rather than pooling it.

Based on the relationship between slope and flow accumulation:
`drainage = sin(slope) / log(1 + accumulation)`
"""
function drainage_potential(dem; method = DInf(), cellsize = cellsize(dem))
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. sind(s) / log1p(acc)
end
