# Single-receiver methods use the same LDD encoding as `flowaccumulation`; methods that can
# send flow several ways need the bitwise D8D encoding.
convention(::D8) = LDD
convention(::DInf) = D8D
convention(::FD8) = D8D

# A method that can send flow several ways needs the share each receiver takes; a single
# receiver takes everything, so there is nothing to carry.
hasweights(::FlowDirectionMethod) = false
hasweights(::DInf) = true
hasweights(::FD8) = true

"""Most receivers a method can ever give one cell."""
maxreceivers(::D8) = 1
maxreceivers(::DInf) = 2
maxreceivers(::FD8) = 8

"""
Share of each cell's flow per receiver, one dimension larger than the DEM, the trailing
axis in [`D8D`](@ref) bit order. A single-receiver method stores nothing: its weights are a
[`Fill`](@extref FillArrays.Fill) and cost no memory.
"""
function allocweights(method::FlowDirectionMethod, dem, T)
    n = maxreceivers(method)
    # A single receiver takes everything, so there is nothing to store.
    hasweights(method) || return Fill(one(T), size(dem)..., n)
    return fill!(_alloc_directions(dem, T, n), zero(T))
end


"""Record the share of a cell's `slot`-th receiver. Nothing to record without weights."""
@inline setshare!(::Fill, cell, slot, fraction) = fraction
@inline setshare!(weights::AbstractArray{T}, cell, slot, fraction) where {T} =
    @inbounds weights[cell, slot] = T(fraction)

"""
Whatever a method wants to set up once for the whole DEM rather than per cell: a grid of
aspects for [`DInf`](@ref) off a raster, and for [`FD8`](@ref) the scratch buffers it
reuses from cell to cell. Methods that need neither get `nothing`.
"""
directionstate(::FlowDirectionMethod, dem, cellsize) = nothing
directionstate(::DInf, dem, cellsize) = _dinf_aspects(dem, cellsize)
directionstate(::FD8, dem, cellsize) = _fd8scratch(dem)

"""
    flowdirection(dem; method=D8(), cellsize=cellsize(dem), weighttype=Float32)

Derive flow directions from `dem` by comparing each cell with its neighbours, and return
them with the share each receiver takes. The encoding matches [`flowaccumulation`](@ref):
[`LDD`](@ref) for the single receiver of [`D8`](@ref), [`D8D`](@ref) for the methods that
can send flow several ways at once.

These depend only on the local neighbourhood, so a cell with no lower neighbour is a pit.
Pass a surface from [`filldepressions`](@ref) if every cell must drain.

# Arguments
- `dem`: An indexed digital elevation model.
- `method::FlowDirectionMethod=D8()`: [`D8`](@ref), [`DInf`](@ref) or [`FD8`](@ref).
- `cellsize`: Cell size, used to scale horizontal distances. Derived from `dem` by default.
- `weighttype=Float32`: Element type of the weights. `Float16` halves the memory a
  divergent method needs, at the cost of a little of the flow, as a cell's shares then no
  longer sum to exactly one. Integer types are not usable.
"""
function flowdirection(
    dem;
    method::FlowDirectionMethod = D8(),
    cellsize = cellsize(dem),
    weighttype::Type{<:AbstractFloat} = Float32,
)
    C = convention(method)
    directions = _directiongrid(dem, C)
    weights = allocweights(method, dem, weighttype)
    state = directionstate(method, dem, cellsize)
    @inbounds for cell in eachcell(dem)
        code = celldirection!(weights, method, state, dem, cell, cellsize)
        directions[cell] = FlowDirection{C}(code)
    end
    return directions, weights
end


"""
Offset towards `cell`'s steepest downslope neighbour, or zero when it has none: the
[`D8`](@ref) rule of
[Jenson and Domingue (1988)](@cite jensonExtractingTopographicStructure1988), and the
fallback wherever a method leaves a cell undecided.
"""
function steepestdescent(dem, cell, cellsize)
    offset = cell - cell
    z = dem[cell]
    isfinite(z) || return offset
    steepest = zero(float(typeof(z)))
    @inbounds for neighbor in neighbors(dem, cell)
        nz = dem[neighbor]
        (isfinite(nz) && nz < z) || continue
        distance = celldistance(dem, cell, neighbor; cellsize)
        iszero(distance) && continue
        gradient = (z - nz) / distance
        if gradient > steepest
            steepest = gradient
            offset = neighbor - cell
        end
    end
    return offset
end

steepestcode(dem, cell, cellsize) =
    _directioncode(dem, steepestdescent(dem, cell, cellsize), cellsize)

celldirection!(weights, ::D8, _, dem, cell, cellsize) =
    _orient(steepestdescent(dem, cell, cellsize), cellsize)

# DInf: the two neighbours bounding the steepest facet, using the same `_dinfreceivers`
# partition as `flowaccumulation`.
function celldirection!(weights, ::DInf, state, dem, cell, cellsize)
    z = dem[cell]
    isfinite(z) || return zero(UInt8)

    receivers = _dinfreceivers(dem, cell, state, cellsize)
    isnothing(receivers) && return steepestshare!(weights, dem, cell, cellsize)
    acell, bcell, afraction, bfraction = receivers

    # The partition can name a cell that is not actually lower, most often beside a flat.
    # Such a receiver takes nothing and the other is renormalised to the whole; if that
    # leaves nobody, fall back rather than strand a cell that can drain.
    lower(ncell) = isfinite(dem[ncell]) && dem[ncell] < z
    a = (_directioncode(dem, acell - cell, cellsize), lower(acell) ? afraction : 0.0)
    b = (_directioncode(dem, bcell - cell, cellsize), lower(bcell) ? bfraction : 0.0)
    total = last(a) + last(b)
    total > 0 || return steepestshare!(weights, dem, cell, cellsize)

    # Shares are read back in the order `decompose` hands the directions out, which is by
    # bit, so emit the pair in that order.
    code = zero(UInt8)
    slot = 0
    @inbounds for (bit, fraction) in (first(a) <= first(b) ? (a, b) : (b, a))
        iszero(fraction) && continue
        slot += 1
        setshare!(weights, cell, slot, fraction / total)
        code |= bit
    end
    return code
end

# FD8: every downslope neighbour receives flow, weighted by slope and contour length, using
# the same `_fd8weights!` as `flowaccumulation`.
function celldirection!(weights, fd8::FD8, scratch, dem, cell, cellsize)
    z = dem[cell]
    isfinite(z) || return zero(UInt8)
    total = _fd8weights!(scratch, fd8, dem, cell, Unvisited(), cellsize)
    iszero(total) && return steepestshare!(weights, dem, cell, cellsize)

    # Gather the shares against their direction bits, then emit them in bit order, which is
    # the order the directions come back out in.
    fill!(scratch.bitweights, 0.0)
    code = zero(UInt8)
    @inbounds for (index, neighbor) in enumerate(scratch.neighbors)
        weight = scratch.weights[index]
        iszero(weight) && continue
        bit = _directioncode(dem, neighbor - cell, cellsize)
        scratch.bitweights[trailing_zeros(Int(bit)) + 1] = weight / total
        code |= bit
    end
    slot = 0
    @inbounds for bit in decompose(D8D, code)
        iszero(bit) && continue
        slot += 1
        setshare!(weights, cell, slot, scratch.bitweights[trailing_zeros(Int(bit)) + 1])
    end
    return code
end

"""The steepest descent fallback: one receiver, taking everything."""
function steepestshare!(weights, dem, cell, cellsize)
    code = steepestcode(dem, cell, cellsize)
    iszero(code) || setshare!(weights, cell, 1, one(Float32))
    return code
end

"""A stencil has no traversal order, so nothing is excluded for being already visited."""
struct Unvisited end
Base.getindex(::Unvisited, _...) = false
