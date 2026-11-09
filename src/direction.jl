# Single-receiver methods use the same LDD encoding as `flowaccumulation`; methods that can
# send flow several ways need the bitwise D8D encoding.
_convention(::D8) = LDD
_convention(::DInf) = D8D
_convention(::FD8) = D8D

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
_weights(method::FlowDirectionMethod, dem, T) =
    hasweights(method) ? zeros(T, size(dem)..., maxreceivers(method)) :
    Fill(one(T), size(dem)..., maxreceivers(method))


"""Record the share of a cell's `slot`-th receiver. Nothing to record without weights."""
@inline _setshare!(::Fill, cell, slot, fraction) = fraction
@inline _setshare!(weights::AbstractArray{T}, cell, slot, fraction) where {T} =
    @inbounds weights[cell, slot] = T(fraction)

# Per-method state cheaper to derive once for the whole DEM than per cell.
_directionstate(::FlowDirectionMethod, dem, cellsize) = nothing
_directionstate(::DInf, dem, cellsize) = _dinf_aspects(dem, cellsize)
_directionstate(::FD8, dem, cellsize) = _fd8scratch(dem)

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
    C = _convention(method)
    directions = similar(dem, FlowDirection{C, UInt8})
    weights = _weights(method, dem, weighttype)
    state = _directionstate(method, dem, cellsize)
    @inbounds for cell in eachcell(dem)
        code = _flowdirection!(weights, method, state, dem, cell, cellsize)
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
function _steepestdescent(dem, cell, cellsize)
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

_steepestcode(dem, cell, cellsize) =
    _directioncode(dem, _steepestdescent(dem, cell, cellsize), cellsize)

_flowdirection!(weights, ::D8, _, dem, cell, cellsize) =
    _orient(_steepestdescent(dem, cell, cellsize), cellsize)

# DInf: the two neighbours bounding the steepest facet, using the same `_dinfreceivers`
# partition as `flowaccumulation`.
function _flowdirection!(weights, ::DInf, state, dem, cell, cellsize)
    z = dem[cell]
    isfinite(z) || return zero(UInt8)

    flowneighbors = _dinfreceivers(dem, cell, state, cellsize)
    isnothing(flowneighbors) && return _steepestshare!(weights, dem, cell, cellsize)
    acell, bcell, afraction, bfraction = flowneighbors

    # A receiver the partition names but that is not lower takes nothing, so the shares of
    # those that remain are renormalised before being recorded.
    takesa = afraction > 0 && isfinite(dem[acell]) && dem[acell] < z
    takesb = bfraction > 0 && isfinite(dem[bcell]) && dem[bcell] < z
    total = (takesa ? afraction : 0.0) + (takesb ? bfraction : 0.0)
    # The partition can name terrain that is not lower, most often beside a flat. Fall back
    # rather than leave a cell that can drain without a direction.
    total > 0 || return _steepestshare!(weights, dem, cell, cellsize)

    abit = takesa ? _directioncode(dem, acell - cell, cellsize) : zero(UInt8)
    bbit = takesb ? _directioncode(dem, bcell - cell, cellsize) : zero(UInt8)
    # Shares are stored in the order the directions come back out, which is by bit
    if iszero(bbit)
        _setshare!(weights, cell, 1, afraction / total)
    elseif iszero(abit)
        _setshare!(weights, cell, 1, bfraction / total)
    elseif abit < bbit
        _setshare!(weights, cell, 1, afraction / total)
        _setshare!(weights, cell, 2, bfraction / total)
    else
        _setshare!(weights, cell, 1, bfraction / total)
        _setshare!(weights, cell, 2, afraction / total)
    end
    return abit | bbit
end

# FD8: every downslope neighbour receives flow, weighted by slope and contour length, using
# the same `_fd8weights!` as `flowaccumulation`.
function _flowdirection!(weights, fd8::FD8, scratch, dem, cell, cellsize)
    z = dem[cell]
    isfinite(z) || return zero(UInt8)
    total = _fd8weights!(scratch, fd8, dem, cell, _Unvisited(), cellsize)
    iszero(total) && return _steepestshare!(weights, dem, cell, cellsize)

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
        _setshare!(weights, cell, slot, scratch.bitweights[trailing_zeros(Int(bit)) + 1])
    end
    return code
end

"""The steepest descent fallback: one receiver, taking everything."""
function _steepestshare!(weights, dem, cell, cellsize)
    code = _steepestcode(dem, cell, cellsize)
    iszero(code) || _setshare!(weights, cell, 1, one(Float32))
    return code
end

"""A stencil has no traversal order, so nothing is excluded for being already visited."""
struct _Unvisited end
Base.getindex(::_Unvisited, _...) = false
