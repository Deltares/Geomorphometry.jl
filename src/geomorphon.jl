"""
    Landform

A landform class of the geomorphon classification, as returned by [`geomorphon`](@ref).
It is stored as an integer, with the same category numbers as GRASS `r.geomorphon`:

- [`Flat`](@ref Landform) (1): no relief in any direction
- [`Peak`](@ref Landform) (2): higher than its surroundings in all directions
- [`Ridge`](@ref Landform) (3): a convex, elongated crest
- [`Shoulder`](@ref Landform) (4): convex break in slope, from flat to falling
- [`Spur`](@ref Landform) (5): convex (nose) slope
- [`Slope`](@ref Landform) (6): a planar, tilted surface
- [`Hollow`](@ref Landform) (7): concave (colluvial) slope
- [`Footslope`](@ref Landform) (8): concave break in slope, from falling to flat
- [`Valley`](@ref Landform) (9): a concave, elongated trough
- [`Pit`](@ref Landform) (10): lower than its surroundings in all directions

[`Undefined`](@ref Landform) (0) marks cells that could not be classified, because their
own elevation is missing or no cell along their lines of sight has a valid elevation.
"""
struct Landform <: Integer
    value::UInt8
    # Inner constructor, so that no default converting constructor is generated that
    # would be ambiguous with the `Integer` constructors in Base.
    Landform(value::Integer) = new(value)
end

const Undefined = Landform(0)
const Flat = Landform(1)
const Peak = Landform(2)
const Ridge = Landform(3)
const Shoulder = Landform(4)
const Spur = Landform(5)
const Slope = Landform(6)
const Hollow = Landform(7)
const Footslope = Landform(8)
const Valley = Landform(9)
const Pit = Landform(10)

const _LANDFORM_NAMES = (
    "Undefined",
    "Flat",
    "Peak",
    "Ridge",
    "Shoulder",
    "Spur",
    "Slope",
    "Hollow",
    "Footslope",
    "Valley",
    "Pit",
)

Base.show(io::IO, l::Landform) = print(io, _LANDFORM_NAMES[l.value + 1])
Base.:(==)(a::Landform, b::Landform) = a.value == b.value
Base.:(<)(a::Landform, b::Landform) = a.value < b.value
Base.Int(l::Landform) = Int(l.value)
Base.UInt8(l::Landform) = l.value
# The category numbers run from `Undefined` to `Pit`; wrapper types use these to pick a
# missing value for a classification.
Base.typemin(::Type{Landform}) = Undefined
Base.typemax(::Type{Landform}) = Pit
# Compare and hash as the category number it stands for.
Base.promote_rule(::Type{Landform}, ::Type{T}) where {T <: Integer} = T
Base.hash(l::Landform, h::UInt) = hash(l.value, h)

# Landform for each (number of negative, number of positive) ternary pattern element
# count, indexed as `_LANDFORMS[nminus + 1, nplus + 1]`. Combinations that sum to more
# than the eight directions are impossible and marked `Undefined`. This is the lookup
# table of [jasiewiczGeomorphonsPatternRecognition2013](@cite).
const _LANDFORMS = @SMatrix [
#   +0        +1        +2        +3        +4        +5        +6        +7        +8
    Flat      Flat      Flat      Footslope Footslope Valley    Valley    Valley    Pit        # -0
    Flat      Flat      Footslope Footslope Footslope Valley    Valley    Valley    Undefined  # -1
    Flat      Shoulder  Slope     Slope     Hollow    Hollow    Valley    Undefined Undefined  # -2
    Shoulder  Shoulder  Slope     Slope     Slope     Hollow    Undefined Undefined Undefined  # -3
    Shoulder  Shoulder  Spur      Slope     Slope     Undefined Undefined Undefined Undefined  # -4
    Ridge     Ridge     Spur      Spur      Undefined Undefined Undefined Undefined Undefined  # -5
    Ridge     Ridge     Ridge     Undefined Undefined Undefined Undefined Undefined Undefined  # -6
    Ridge     Ridge     Undefined Undefined Undefined Undefined Undefined Undefined Undefined  # -7
    Peak      Undefined Undefined Undefined Undefined Undefined Undefined Undefined Undefined  # -8
]

"""
    geomorphon(dem; cellsize=cellsize(dem), radius=10, skip=0, flatness=1.0)

Classify each cell of a DEM into one of ten [`Landform`](@ref) classes, using the
geomorphon approach of [Jasiewicz and Stepinski (2013)](@cite jasiewiczGeomorphonsPatternRecognition2013).

Along each of the eight principal directions the highest (zenith) and lowest (nadir)
elevation angles to the cells within the search radius are determined. Their sum is
positive where the terrain in that direction rises above the cell, and negative where it
falls below, which reduces each direction to `+`, `0` or `-`. The resulting ternary
pattern of eight elements is summarized by the number of `+` and `-` elements, which maps
onto the ten most common landforms.

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `radius`: Outer search radius in cells (default `10`), measured from cell centre to cell
  centre, so the sampled cells are those of a `Circle(radius)` neighborhood.
- `skip`: Inner search radius in cells (default `0`). The first `skip` cells of each line
  of sight are ignored, to suppress local noise.
- `flatness`: Flatness threshold in degrees (default `1.0`). Directions whose zenith and
  nadir angles sum to less than this are considered flat.

# Returns
A matrix of [`Landform`](@ref) values. Cells with a `NaN` elevation are `Undefined`; cells
with missing elevations along a line of sight are skipped.
"""
function geomorphon(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    radius::Integer = 10,
    skip::Integer = 0,
    flatness::Real = 1.0,
)
    radius >= 1 || throw(ArgumentError("radius must be at least 1, got $radius"))
    0 <= skip < radius || throw(ArgumentError("skip must be in 0:$(radius - 1), got $skip"))
    flatness >= 0 || throw(ArgumentError("flatness must not be negative, got $flatness"))

    cs1, cs2 = abs(Float64(cellsize[1])), abs(Float64(cellsize[2]))
    nrows, ncols = size(dem)

    out = _alloc(dem, Landform)
    backend = get_backend(parent(dem))

    workgroup = backend isa KernelAbstractions.CPU ? 1 : (16, 16)
    kernel! = _geomorphon_kernel!(backend, workgroup)
    kernel!(
        out,
        dem,
        cs1,
        cs2,
        Int(radius),
        Int(skip),
        Float64(flatness),
        nrows,
        ncols;
        ndrange = size(dem),
    )
    KernelAbstractions.synchronize(backend)

    return out
end

@kernel function _geomorphon_kernel!(
    out,
    @Const(dem),
    cs1::Float64,
    cs2::Float64,
    radius::Int,
    skip::Int,
    flatness::Float64,
    nrows::Int,
    ncols::Int,
)
    i, j = @index(Global, NTuple)
    out[i, j] =
        _geomorphon(dem, i, j, cs1, cs2, radius, skip, flatness, nrows, ncols)
end

# Reduce the eight lines of sight around (i, j) to a ternary pattern and look up its
# landform. Only the counts of rising (+) and falling (-) directions matter, as the
# lookup is invariant under rotation and mirroring of the pattern.
@inline function _geomorphon(
    dem,
    i,
    j,
    cs1,
    cs2,
    radius,
    skip,
    flatness,
    nrows,
    ncols,
)
    @inbounds z = Float64(dem[i, j])
    isnan(z) && return Undefined

    nplus = 0
    nminus = 0
    @inbounds for dr in -1:1, dc in -1:1
        (dr == 0 && dc == 0) && continue
        # One step along this ray, in cells and in ground units; the search radius counts
        # cells, while the angles are measured over the real distance.
        step = hypot(dr, dc)
        metric_step = hypot(dr * cs1, dc * cs2)
        # Highest (zenith) and lowest (nadir) elevation angle along the ray, kept as
        # tangents; the tangent is monotonic in the angle, so the extremes can be
        # tracked without converting each sample.
        zenith = -Inf
        nadir = Inf
        k = skip + 1
        while k * step <= radius
            r = i + k * dr
            c = j + k * dc
            (r < 1 || r > nrows || c < 1 || c > ncols) && break
            distance = k * metric_step
            zk = Float64(dem[r, c])
            if !isnan(zk)
                tangent = (zk - z) / distance
                zenith = max(zenith, tangent)
                nadir = min(nadir, tangent)
            end
            k += 1
        end
        zenith == -Inf && continue  # no valid cell in this direction

        # The zenith angle is measured up and the nadir angle down from the cell, so
        # their sum is the excess of terrain above it.
        excess = atand(zenith) + atand(nadir)
        if excess > flatness
            nplus += 1
        elseif excess < -flatness
            nminus += 1
        end
    end

    return @inbounds _LANDFORMS[nminus + 1, nplus + 1]
end
