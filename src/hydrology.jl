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

"""Priority Flood."""
function priorityflood(dem::AbstractMatrix, mask = falses(size(dem)))
    priorityflood!(copy(dem), mask)
end

const Δ = CartesianIndex(1, 1)

abstract type FlowDirectionMethod end

"""D8 Flow Direction method (Jenson and Domingue, 1988)."""
struct D8 <: FlowDirectionMethod end

"""DInf Flow Direction method (Tarboton, 1997)."""
struct DInf <: FlowDirectionMethod end

function priorityflood!(dem::AbstractMatrix, closed = falses(size(dem)))
    open = PriorityQueue{CartesianIndex{2}, eltype(dem)}()
    pit = DataStructures.Queue{CartesianIndex{2}}()

    R = CartesianIndices(dem)
    I_first, I_last = first(R), last(R)

    @inbounds for cell in edges(R)
        enqueue!(open, cell, dem[cell])
        closed[cell] = true
    end
    @inbounds while !isempty(open) || !isempty(pit)
        cell = !isempty(pit) ? DataStructures.dequeue!(pit) : dequeue!(open)
        for ncell in max(I_first, cell - Δ):min(I_last, cell + Δ)
            (closed[ncell] || ncell == cell) && continue
            closed[ncell] = true
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

"""Streamflow."""
function streamflow(dem::AbstractMatrix, closed = falses(size(dem)); method = D8())
    streamflow!(dem, copy(closed); method)
end

function streamflow!(dem::AbstractMatrix, closed = falses(size(dem)); method)
    asp = aspect(dem; method = Horn())

    dir = fill(CartesianIndex{2}(0, 0), size(dem))
    order = ones(Int64, length(closed) - sum(closed))

    acc = similar(dem, Float32)
    acc .= 1
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

    _accumulate!(method, acc, asp, order, dir, R)
    output .= getindex.(Ref(directions), dir)
    return acc, output
end

function _accumulate!(::D8, acc, asp, order, dir, R)
    for i in reverse(order)
        acc[R[i] + dir[i]] += acc[i]
    end
end
function _accumulate!(::DInf, acc, asp, order, dir, R)
    for i in reverse(order)
        aspect = asp[i]
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
