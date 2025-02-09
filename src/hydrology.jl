using OffsetArrays: centered

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

Performs the Priority-Flood algorithm on the given digital elevation model (DEM) `dem` with an optional `mask`.

# Arguments
- `dem::AbstractMatrix`: A 2D array representing the digital elevation model (DEM).
- `mask::AbstractMatrix{Bool}`: A 2D boolean array representing the mask. Cells with `true` values are considered in the computation, while cells with `false` values are ignored.
"""
function filldepressions(dem::AbstractMatrix, mask = falses(size(dem)))
    filldepressions!(copy(dem), mask)
end

const Δ = CartesianIndex(1, 1)

abstract type FlowDirectionMethod end

"""D8 Flow Direction method (Jenson and Domingue, 1988)."""
struct D8 <: FlowDirectionMethod end

"""DInf Flow Direction method (Tarboton, 1997)."""
struct DInf <: FlowDirectionMethod end

"""FD8 Flow Direction method (Quinn et al., 1991)."""
Base.@kwdef struct FD8 <: FlowDirectionMethod
    L1::Float32 = 0.5
    L2::Float32 = 0.354
    p::Float32 = 1
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
        if queued[cell] && !ismissing(dem[cell]) && isfinite(dem[cell]) && iszero(labels[cell])
            labels[cell] = label
            label += 1
        end
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
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

const nb =
    CartesianIndex.([
        (0, -1),
        (-1, -1),
        (-1, 0),
        (-1, 1),
        (0, 1),
        (1, 1),
        (1, 0),
        (1, -1),
        (0, -1),
        (-1, -1),
    ])

function infc(aspect)
    i = round(Int, aspect / 45) + 1
    return nb[i], nb[i + 1]
end

function infa(aspect)
    a1 = aspect % 45
    a2 = 45 - a1
    return a2 / 45, a1 / 45
end

const directions = centered([1 2 3; 4 5 6; 7 8 9])

"""
    flowaccumulation(dem::AbstractMatrix, closed::Matrix{Bool}, method::FlowDirectionMethod)

Computes the flow accumulation of a digital elevation model (DEM) `dem` with an optional `closed` mask and a `method` for flow direction.
Returns the flow accumulation and the flow direction (local drainage direction or ldd)
"""
function flowaccumulation(dem::AbstractMatrix, closed = falses(size(dem)); method = D8(), cellsize=1)
    flowaccumulation!(dem, copy(closed); method, cellsize)
end

function flowaccumulation!(dem::AbstractMatrix, closed = falses(size(dem)); method=D8(), cellsize=1)

    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    order = ones(Int64, length(closed) - sum(closed))

    acc = similar(dem, Float32)
    acc .= cellsize^2
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
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            # skip visited and center cells
            (closed[ncell] || ncell == cell) && continue

            closed[ncell] = true
            dir[ncell] = cell - ncell

            enqueue!(open, ncell, dem[ncell])
        end
    end

    _accumulate!(method, acc, order, dir, R, dem)
    output .= getindex.(Ref(directions), dir)
    return acc, output
end

function _accumulate!(::D8, acc, order, dir, R, dem)
    for i in reverse(order)
        acc[R[i] + dir[i]] += acc[i]
    end
end
function _accumulate!(::DInf, acc, order, dir, R, dem)
    asp = aspect(dem; method = Horn())
    for i in reverse(order)
        aspect = asp[i]

        if !isfinite(aspect)
            acc[R[i] + dir[i]] += acc[i]
            continue
        end

        a, b = infc(aspect)
        aa, ab = infa(aspect)

        # Depression
        if a != dir[i] && b != dir[i]
            acc[R[i] + dir[i]] += acc[i]
            continue
        end

        if R[i] + a in R
            acc[R[i] + a] += acc[i] * aa
        end
        if R[i] + b in R
            acc[R[i] + b] += acc[i] * ab
        end
    end
end

function _accumulate!(fd8::FD8, acc, order, dir, R, dem)
    contour_lengths = @SMatrix [
        fd8.L2 fd8.L1 fd8.L2
        fd8.L1 0 fd8.L1
        fd8.L2 fd8.L1 fd8.L2
    ]
    weights = zeros(size(contour_lengths))
    Σw = 0.0
    for i in reverse(order)
        Σw = 0.0
        for (ri, dist) in enumerate(neib_8_dist)
            ri == 5 && continue
            I = R[i] + nb[ri]
            I in R || continue
            diff = dem[R[i]] - dem[I]
            if diff < 0 || isnan(diff) # neighbor is higher
                weights[ri] = 0.0
                continue
            end
            weight = (diff / dist * contour_lengths[ri])^fd8.p
            weights[ri] = weight
            Σw += weight
        end
        if iszero(Σw)  # depression
            acc[R[i] + dir[i]] += acc[i]
            continue
        end
        for (ri, weight) in enumerate(weights)
            iszero(weight) && continue
            I = R[i] + nb[ri]
            I in R || continue
            acc[I] += acc[i] * (weight / Σw)
        end
    end
end

"""
    TWI(dem::AbstractMatrix; method=D8(), cellsize=1)

Computes the Topographic Wetness Index (TWI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function TWI(dem::AbstractMatrix; method=D8(), cellsize=1)
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc / tand(s))
end

"""
    SPI(dem::AbstractMatrix; method=D8(), cellsize=1)

Computes the Stream Power Index (SPI) of a digital elevation model (DEM) `dem` with an optional `method` for flow direction and a `cellsize`.
"""
function SPI(dem::AbstractMatrix; method=D8(), cellsize=1)
    s = slope(dem; cellsize)
    acc, _ = flowaccumulation(dem; method, cellsize)
    return @. log(acc * tand(s))
end
