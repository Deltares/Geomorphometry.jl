module GeomorphometryDiscreteGlobalGridsExt

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import Rasters
using Rasters: Raster

# One-dimensional rasters indexed by a DGG `Cells` lookup.
const CellsRaster{T,D} =
    Raster{T,1,D} where {T,D<:Tuple{<:DGG.Cells{<:DGG.CellLookup}}}
const IGeo7Raster{T,D} =
    Raster{T,1,D} where {T,D<:Tuple{<:DGG.Cells{<:DGG.CellLookup{DGG.Z7Cell}}}}
const RelativeIndex = DGG.RelativeZ7Cell
const Cell = DGG.Z7Cell
const CellIndex = DGG.AbstractCellIndex

# Authalic Earth radius used to scale unit-sphere measurements.
const R_AUTHALIC = DGG.ISEA.R_AUTHALIC

_lookup(r::CellsRaster) = Rasters.lookup(r, DGG.Cells)
_cellvector(r::CellsRaster) = parent(_lookup(r))
_completegrid(cells) = DGG.levelgrid(DGG.system(cells), DGG.level(cells))

# `CellVector` is lazy and its `in` method searches compressed position windows
# rather than scanning cell ids.
Base.eachindex(r::CellsRaster) = _cellvector(r)

function _position(r::CellsRaster, c::CellIndex)
    p = DGG.cellposition(_lookup(r), c)
    # Wrap the cell because `BoundsError` iterates its index field.
    isnothing(p) && throw(BoundsError(r, (c,)))
    return p
end

Base.getindex(r::CellsRaster, c::CellIndex) = parent(r)[_position(r, c)]
Base.setindex!(r::CellsRaster, value, c::CellIndex) =
    setindex!(parent(r), value, _position(r, c))
Base.checkbounds(::Type{Bool}, r::CellsRaster, c::CellIndex) =
    DGG.cellposition(_lookup(r), c) !== nothing

GM.neighbors(r::CellsRaster, c::CellIndex) = DGG.neighbors(_lookup(r), c)

# Neighborhood traversal

GM.neighbors(r::CellsRaster) = DGG.neighbors(r)

# Let DGG resolve the cell axis and rebuild the raster. `Values()` passes the
# cell and neighbor values to `f`.
GM.mapneighbors(f::F, r::CellsRaster; order = nothing, threaded = true) where {F} =
    DGG.mapneighbors(f, r; pass = DGG.Values(),
        order = order === nothing ? DGG.StorageOrder() : order, threaded)

# Boundary cells have fewer neighbors than on the complete level grid.
function GM.outlets(r::CellsRaster)
    cells = _cellvector(r)
    complete = _completegrid(cells)
    rim = DGG.mapneighbors(
        (c, nbrs) -> length(nbrs) < DGG.neighborcount(complete, DGG.cellid(c)),
        cells; threaded = true)
    return [cells[p] for p in eachindex(rim) if rim[p]]
end

# Priority-flood traversal

# Visit cells from the boundary inward, taking the lowest queued elevation.
# `order` records visited positions and `down[p]` records the position that
# first reached `p`. Closed cells are skipped, except that boundary cells are
# always queued.
function _settle!(order::Vector{Int64}, down::Vector{Int}, closedv, zv, cv,
        table)
    n = length(cv)
    complete = _completegrid(cv)
    open = GM.FastPriorityQueue{eltype(zv)}(n)
    @inbounds for p in 1:n
        if length(table[p]) < DGG.neighborcount(complete, cv[p])
            GM.enqueue!(open, p, zv[p])
            closedv[p] = true
        end
    end
    i = 1
    @inbounds while !isempty(open)
        p = GM.dequeue!(open)
        order[i] = p
        i += 1
        for q in table[p]
            closedv[q] && continue
            closedv[q] = true
            down[q] = p
            GM.enqueue!(open, q, zv[q])
        end
    end
    return nothing
end

# Convert downstream positions to relative IGeo7 directions; zero marks no outflow.
function _directions(dem::IGeo7Raster, cv, down)
    zrel = first(cv) - first(cv)
    dir = similar(dem, typeof(zrel); missingval = nothing)
    dirv = parent(dir)
    @inbounds for p in eachindex(down)
        dirv[p] = down[p] == 0 ? zrel : cv[down[p]] - cv[p]
    end
    return dir
end

# Match the missingval convention of the D8 direction outputs below; the
# direction eltype cannot represent the DEM missingval anyway.
GM._directiongrid(dem::CellsRaster, ::Type{C}) where {C <: GM.FlowDirectionConvention} =
    similar(dem, GM.FlowDirection{C, UInt8}; missingval = nothing)

# IGeo7 uses its constant cell area; other grids use their spherical polygon area.
_cellarea(grid, c::CellIndex) = DGG.cell_area(grid, c) * R_AUTHALIC^2
_cellarea(grid, c::Cell) = DGG.IGeo7.cell_area(DGG.rawid(c))

function GM.flowaccumulation(dem::CellsRaster,
        closed = fill!(similar(dem, Bool; missingval = nothing), false);
        method = GM.D8(), cellsize = GM.cellsize(dem))
    cv = _cellvector(dem)
    grid = _completegrid(cv)
    acc = similar(dem, Float32; missingval = NaN32)
    accv = parent(acc)
    @inbounds for p in eachindex(accv)
        accv[p] = _cellarea(grid, cv[p])
    end
    return GM.flowaccumulation!(dem, acc, copy(closed); method, cellsize)
end

function GM.flowaccumulation!(dem::CellsRaster, acc::AbstractArray{<:Real},
        closed = fill!(similar(dem, Bool; missingval = nothing), false);
        method = GM.D8(), cellsize = GM.cellsize(dem))
    cv = _cellvector(dem)
    closedv = parent(closed)
    order = ones(Int64, length(cv) - count(closedv))
    down = zeros(Int, length(cv))
    table = DGG.adjacency(cv)
    _settle!(order, down, closedv, parent(dem), cv, table)
    return _postsettle(method, dem, acc, order, down, cv, cellsize)
end

# Accumulate each cell into its downstream position in reverse traversal order.
function _accumulate_down!(accv, order, down)
    @inbounds for j in reverse(eachindex(order))
        p = order[j]
        d = down[p]
        d == 0 && continue
        accv[d] += accv[p]
    end
    return accv
end

# Flow direction encoding
#
# Non-IGeo7 grids encode the downstream neighbor's complete-level ring slot as
# a `UInt16` bit. Slot `k` sets bit `k - 1`, and zero marks no outflow. Vertex
# rings include corner neighbors on quad-based grids and may vary in length.

function _ringslot(complete, cell::CellIndex, target::CellIndex)
    ring = DGG.neighbors(complete, cell)
    slot = findfirst(==(target), ring)
    slot === nothing &&
        throw(ArgumentError("$target is not a ring member of $cell"))
    slot <= 16 ||
        throw(ArgumentError("ring slot $slot does not fit the UInt16 codec"))
    return UInt16(1) << (slot - 1)
end

# Decode set slots to cells in ring order.
function _slottargets(complete, cell::CellIndex, d::GM.FlowDirection{GM.D8D})
    ring = DGG.neighbors(complete, cell)
    bits = Int(d)
    return (ring[k] for k in 1:length(ring) if !iszero(bits & (1 << (k - 1))))
end

function _postsettle(::GM.D8, dem::CellsRaster, acc, order, down, cv, cellsize)
    _accumulate_down!(parent(acc), order, down)
    complete = _completegrid(cv)
    output = similar(dem, GM.FlowDirection{GM.D8D, UInt16}; missingval = nothing)
    outv = parent(output)
    @inbounds for p in eachindex(down)
        outv[p] = GM.FlowDirection{GM.D8D}(down[p] == 0 ? UInt16(0) :
            _ringslot(complete, cv[p], cv[down[p]]))
    end
    return acc, output
end

_postsettle(method::GM.FlowDirectionMethod, dem::CellsRaster, acc, order,
    down, cv, cellsize) = throw(ArgumentError(
    "flow accumulation with $(nameof(typeof(method))) needs relative-cell " *
    "arithmetic, which only the IGeo7 backend provides; use D8"))

# IGeo7 stores D8 results with its relative-cell LDD encoding.
function _postsettle(::GM.D8, dem::IGeo7Raster, acc, order, down, cv, cellsize)
    _accumulate_down!(parent(acc), order, down)
    zrel = first(cv) - first(cv)
    output = similar(dem, GM.FlowDirection{GM.LDD, UInt8}; missingval = nothing)
    outv = parent(output)
    @inbounds for p in eachindex(down)
        rel = down[p] == 0 ? zrel : cv[down[p]] - cv[p]
        outv[p] = GM.FlowDirection{GM.LDD}(rel)
    end
    return acc, output
end

# DInf and FD8 derive flow shares from the DEM and relative IGeo7 directions.
function _postsettle(method::GM.FlowDirectionMethod, dem::IGeo7Raster, acc,
        order, down, cv, cellsize)
    dir = _directions(dem, cv, down)
    dirs = GM._accumulate!(method, acc, order, dir, cv, dem, cellsize)
    return acc, dirs
end

function GM.height_above_nearest_drainage(dem::IGeo7Raster;
        method = GM.D8(), cellsize = GM.cellsize(dem), threshold = 100)
    cv = _cellvector(dem)
    n = length(cv)
    output = zero(dem)
    acc = similar(dem, Float32; missingval = NaN32)
    accv = parent(acc)
    @inbounds for p in 1:n
        accv[p] = DGG.IGeo7.cell_area(DGG.rawid(cv[p]))
    end
    closed = fill!(similar(dem, Bool; missingval = nothing), false)
    order = ones(Int64, n)
    down = zeros(Int, n)
    table = DGG.adjacency(cv)
    _settle!(order, down, parent(closed), parent(dem), cv, table)
    dir = _directions(dem, cv, down)
    flowdirs = GM._accumulate!(method, acc, order, dir, cv, dem, cellsize)
    stream_mask = acc .>= threshold
    GM._burn_streams!(stream_mask, order, dir, cv)
    GM._fill_flow_depressions!(acc, order, flowdirs, cv, dem, cellsize)
    GM._hand!(output, order, flowdirs, cv, acc, stream_mask, cellsize)
    return output
end
      
function GM.outlets(r::IGeo7Raster)
    cells = _cells(r)
    complete = DGG.levelgrid(DGG.system(cells), DGG.level(cells))
    return Cell[
        c for c in cells
              if length(DGG.neighbors(cells, c)) < length(DGG.neighbors(complete, c))
    ]
end

# DGGS HAND supports only the threshold form on IGeo7 rasters. Reject the
# other entry points before they reach position-based code.
GM.height_above_nearest_drainage(dem::CellsRaster; kw...) =
    throw(ArgumentError(
        "height_above_nearest_drainage needs relative-cell arithmetic, " *
        "which only the IGeo7 backend provides"))

GM.height_above_nearest_drainage(dem::CellsRaster,
    stream_mask::AbstractArray{Bool}) = throw(ArgumentError(
    "height_above_nearest_drainage with a stream mask is not implemented " *
    "for DGGS rasters; use the threshold form, which the IGeo7 backend " *
    "provides"))

function GM.cellarea(r::CellsRaster, c::CellIndex; cellsize=nothing)
    _position(r, c)
    return _cellarea(DGG.levelgrid(DGG.system(_cellvector(r)), DGG.level(c)), c)
end

function GM.celldistance(
    r::CellsRaster,
    from::CellIndex,
    to::CellIndex;
    cellsize=nothing,
)
    _position(r, from)
    _position(r, to)
    grid = DGG.levelgrid(DGG.system(_cellvector(r)), DGG.level(from))
    a = DGG.cell_centroid(grid, from)
    b = DGG.cell_centroid(grid, to)
    angle = acos(clamp(a[1] * b[1] + a[2] * b[2] + a[3] * b[3], -1.0, 1.0))
    return angle * R_AUTHALIC
end

function GM.cellbearing(r::IGeo7Raster, from::Cell, to::Cell; cellsize=nothing)
    _position(r, from)
    _position(r, to)
    from == to && return 0.0
    grid = DGG.levelgrid(DGG.system(_cellvector(r)), DGG.level(from))
    a = DGG.cell_centroid(grid, from)
    b = DGG.cell_centroid(grid, to)
    lon1, lon2 = atan(a[2], a[1]), atan(b[2], b[1])
    lat1, lat2 = asin(clamp(a[3], -1.0, 1.0)), asin(clamp(b[3], -1.0, 1.0))
    delta_lon = lon2 - lon1
    east = sin(delta_lon) * cos(lat2)
    north = cos(lat1) * sin(lat2) -
            sin(lat1) * cos(lat2) * cos(delta_lon)
    return mod(rad2deg(atan(east, north)), 360.0)
end

const IGEO7_TO_LDD = (0x05, 0x06, 0x09, 0x07, 0x04, 0x01, 0x03)
const IGEO7_TO_D8D = (0x00, 0x01, 0x80, 0x20, 0x10, 0x08, 0x02)

GM.FlowDirection{GM.LDD}(d::RelativeIndex) =
    GM.FlowDirection{GM.LDD}(IGEO7_TO_LDD[DGG.directioncode(d)+1])
GM.FlowDirection{GM.D8D}(d::RelativeIndex) =
    GM.FlowDirection{GM.D8D}(IGEO7_TO_D8D[DGG.directioncode(d)+1])

function _inverse_code(codes, value)
    index = findfirst(==(UInt8(value)), codes)
    return isnothing(index) ? 0xff : UInt8(index - 1)
end

const LDD_TO_IGEO7 = ntuple(value -> _inverse_code(IGEO7_TO_LDD, value), 9)
const D8D_BIT_TO_IGEO7 =
    ntuple(bit -> _inverse_code(IGEO7_TO_D8D, 1 << (bit - 1)), 8)

@inline function _directioncode(direction::GM.FlowDirection{GM.LDD})
    value = Int(direction)
    1 <= value <= length(LDD_TO_IGEO7) ||
        throw(ArgumentError("invalid LDD direction $value"))
    code = @inbounds LDD_TO_IGEO7[value]
    code != 0xff ||
        throw(ArgumentError("LDD direction $value has no IGeo7 equivalent"))
    return code
end

@inline function _directioncode(direction::GM.FlowDirection{GM.D8D})
    iszero(direction) && return 0x00
    value = Int(direction)
    ispow2(value) ||
        throw(ArgumentError("D8D direction must be decomposed before conversion"))
    bit = trailing_zeros(value) + 1
    bit <= length(D8D_BIT_TO_IGEO7) ||
        throw(ArgumentError("D8D direction $value has no IGeo7 equivalent"))
    code = @inbounds D8D_BIT_TO_IGEO7[bit]
    code != 0xff ||
        throw(ArgumentError("D8D direction $value has no IGeo7 equivalent"))
    return code
end

function GM.decompose(
    ::Type{RelativeIndex},
    direction::GM.FlowDirection,
    center::Cell,
)
    return Iterators.map(GM.decompose(direction)) do component
        code = _directioncode(component)
        iszero(code) && return RelativeIndex(center, 0)
        adjacent = DGG.neighbors(
            DGG.levelgrid(DGG.IGeo7System(), DGG.level(center)),
            center,
        )
        code <= length(adjacent) ||
            throw(ArgumentError("direction $code does not exist at pentagon $center"))
        return adjacent[code] - center
    end
end

end # module
