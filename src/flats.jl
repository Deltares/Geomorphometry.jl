"""
    resolveflats(dem)
    resolveflats!(dem)

Give the flats of `dem` a drainage gradient, following
[Barnes et al. (2014)](@cite barnesEfficientAssignmentDrainage2014), and return the
elevation model. `resolveflats!` writes into `dem` instead of a copy.

A flat is a connected region of equal elevation, over which no flow direction is defined.
The gradient is the superposition of one away from higher terrain, which gives convergent
flow, and one towards lower terrain, which guarantees the flat drains. The latter carries
twice the weight, so it settles cells where the two disagree without iterating.

Cells are raised by whole `nextfloat` steps, so the surface is unchanged at any meaningful
precision but every flat cell gains a strictly lower neighbour. A flat with no outlet is a
depression, not a flat, and is left alone. [`filldepressions`](@ref) calls this by default.
"""
resolveflats(dem) = resolveflats!(copy(dem))

function resolveflats!(dem)
    increments = _flatmask(dem)
    @inbounds for cell in eachcell(dem)
        increment = increments[cell]
        iszero(increment) || (dem[cell] = nextfloat(dem[cell], increment))
    end
    return dem
end

"""
Per-cell increment counts for [`resolveflats`](@ref): zero wherever a cell already drains,
positive on cells that have no strictly lower neighbour and belong to a flat that drains.
"""
function _flatmask(dem)
    R = eachcell(dem)
    Cell = typeof(first(R))

    drains = fill!(similar(dem, Bool), false)
    labels = fill!(similar(dem, Int32), 0)
    flatmask = fill!(similar(dem, Int32), 0)

    highedges = Cell[]
    lowedges = Cell[]

    # Cells on the boundary drain out of the grid, so they never need resolving.
    @inbounds for cell in outlets(dem)
        drains[cell] = true
    end

    # Step 1: a cell drains if it has a strictly lower neighbour; a high edge cell cannot
    # drain but touches higher terrain.
    @inbounds for cell in R
        z = dem[cell]
        if !isfinite(z)
            drains[cell] = true
            continue
        end
        drains[cell] && continue
        lower = false
        higher = false
        for ncell in neighbors(dem, cell)
            nz = dem[ncell]
            isfinite(nz) || continue
            nz < z && (lower = true)
            nz > z && (higher = true)
        end
        drains[cell] = lower
        (!lower && higher) && push!(highedges, cell)
    end

    # A low edge cell drains, but borders a cell of its own elevation that does not: it is
    # an outlet of that flat.
    @inbounds for cell in R
        drains[cell] || continue
        isfinite(dem[cell]) || continue
        for ncell in neighbors(dem, cell)
            if dem[ncell] == dem[cell] && !drains[ncell]
                push!(lowedges, cell)
                break
            end
        end
    end
    isempty(lowedges) && return flatmask

    # Label each flat by flood filling its outlets over cells of equal elevation. Flats
    # without an outlet keep label 0 and are left alone: they are depressions.
    nlabels = zero(Int32)
    queue = Cell[]
    @inbounds for seed in lowedges
        iszero(labels[seed]) || continue
        nlabels += one(nlabels)
        z = dem[seed]
        labels[seed] = nlabels
        empty!(queue)
        push!(queue, seed)
        head = 1
        while head <= length(queue)
            cell = queue[head]
            head += 1
            for ncell in neighbors(dem, cell)
                if iszero(labels[ncell]) && dem[ncell] == z
                    labels[ncell] = nlabels
                    push!(queue, ncell)
                end
            end
        end
    end
    filter!(cell -> !iszero(labels[cell]), highedges)

    # Step 2: gradient away from higher terrain, breadth-first inwards from the high
    # edges. `flatheights` is how far each flat reaches, which Step 3 needs to invert it.
    flatheights = zeros(Int32, nlabels)
    seen = fill!(similar(dem, Bool), false)
    _seed!(seen, highedges)
    increment = one(Int32)
    frontier = copy(highedges)
    next = Cell[]
    @inbounds while !isempty(frontier)
        empty!(next)
        for cell in frontier
            flatmask[cell] = increment
            flatheights[labels[cell]] = increment
            _growflat!(next, seen, dem, labels, drains, cell)
        end
        frontier, next = next, frontier
        increment += one(increment)
    end

    # Step 3: gradient towards lower terrain, growing back from the outlets. Negating
    # Step 2 and adding the flat's height inverts it into the gradient we want; every cell
    # then gains twice its distance from the outlet, so draining dominates on disagreement.
    @inbounds for cell in R
        flatmask[cell] = -flatmask[cell]
    end
    fill!(seen, false)
    _seed!(seen, lowedges)
    increment = one(Int32)
    frontier = copy(lowedges)
    empty!(next)
    @inbounds while !isempty(frontier)
        empty!(next)
        for cell in frontier
            _growflat!(next, seen, dem, labels, drains, cell)
            flatmask[cell] < 0 && (flatmask[cell] += flatheights[labels[cell]])
            flatmask[cell] += 2 * increment
        end
        frontier, next = next, frontier
        increment += one(increment)
    end

    # Cells that already drain keep their elevation: the gradient is only for cells that
    # would otherwise have no flow direction.
    @inbounds for cell in R
        drains[cell] && (flatmask[cell] = 0)
    end
    return flatmask
end

function _seed!(seen, cells)
    for cell in cells
        seen[cell] = true
    end
    return seen
end

"""Queue the not-yet-visited cells of `cell`'s own flat that still need a gradient."""
@inline function _growflat!(next, seen, dem, labels, drains, cell)
    label = labels[cell]
    @inbounds for ncell in neighbors(dem, cell)
        if !seen[ncell] && labels[ncell] == label && !drains[ncell]
            seen[ncell] = true
            push!(next, ncell)
        end
    end
    return next
end
