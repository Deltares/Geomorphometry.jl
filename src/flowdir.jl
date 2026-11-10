"""
    FlowDirectionConvention

Abstract type for flow direction encoding conventions. Subtypes define how
direction integers map to neighbor offsets. See [`LDD`](@ref).
"""
abstract type FlowDirectionConvention end

"""
    FlowDirection{C<:FlowDirectionConvention, T<:Integer} <: Integer

A flow direction value in convention `C`, stored as type `T`.

# Constructors
- `FlowDirection{C}(v::Integer)`: Create from a raw direction value.
- `FlowDirection{C}(ci::CartesianIndex{2})`: Create from a CartesianIndex offset.

# Conversion
- `CartesianIndex(d::FlowDirection)`: Convert to a CartesianIndex offset.
- `decompose(R, d::FlowDirection, center)`: Decompose into relative indices of type `R`.
- `convert(FlowDirection{C2}, d::FlowDirection{C1})`: Convert between conventions.
- `Int(d::FlowDirection)`: Get the raw integer value.

Grid implementations can support flow-direction traversal by defining
`decompose(::Type{R}, d::FlowDirection, center)` for their relative-index type `R`.
"""
struct FlowDirection{C <: FlowDirectionConvention, T <: Integer} <: Integer
    value::T
    FlowDirection{C, T}(value::Integer) where {C <: FlowDirectionConvention, T <: Integer} =
        new{C, T}(value)
end

FlowDirection{C}(v::T) where {C <: FlowDirectionConvention, T <: Integer} =
    FlowDirection{C, T}(v)
FlowDirection{C}(ci::CartesianIndex{2}) where {C <: FlowDirectionConvention} =
    FlowDirection{C}(_ci_to_dir(C, ci))

Base.show(io::IO, d::FlowDirection{C}) where {C} = print(io, _arrow(C, d.value))
Base.CartesianIndex(d::FlowDirection{C}) where {C} = _dir_to_ci(C, d.value)
Base.Int(d::FlowDirection) = Int(d.value)
Base.:(==)(a::FlowDirection{C}, b::FlowDirection{C}) where {C} = a.value == b.value
Base.:(==)(a::FlowDirection{C1}, b::FlowDirection{C2}) where {C1, C2} = false
Base.:(==)(a::FlowDirection{C}, b::Integer) where {C} = a.value == b
Base.:(==)(a::Integer, b::FlowDirection{C}) where {C} = a == b.value
Base.:(==)(a::FlowDirection{C}, b::BigInt) where {C} = a.value == b
Base.:(==)(a::BigInt, b::FlowDirection{C}) where {C} = a == b.value

Base.convert(::Type{CartesianIndex{2}}, d::FlowDirection{C}) where {C} =
    _dir_to_ci(C, d.value)
Base.convert(::Type{FlowDirection{C}}, ci::CartesianIndex{2}) where {C} =
    FlowDirection{C}(_ci_to_dir(C, ci))
Base.convert(t::Type{T}, d::FlowDirection{C, T}) where {C, T <: Integer} = d.value
Base.convert(::Type{FlowDirection{C}}, d::T) where {C, T <: Integer} =
    FlowDirection{C, T}(d)
Base.convert(::Type{FlowDirection{C}}, d::FlowDirection{C2}) where {C, C2} =
    FlowDirection{C}(Base.convert(CartesianIndex{2}, d))
Base.convert(::Type{FlowDirection{C}}, d::FlowDirection{C, <:FlowDirection}) where {C} = d

"""Return the `FlowDirectionConvention` type of a direction."""
convention(::Type{FlowDirection{C}}) where {C} = C
convention(::Type{FlowDirection{C, T}}) where {C, T} = C
convention(d::FlowDirection) = convention(typeof(d))

"""Whether a convention supports encoding multiple directions in a single value."""
ismulti(::Type{<:FlowDirectionConvention}) = false

"""
    LDD <: FlowDirectionConvention

Local Drainage Direction (PCRaster) convention using 1-9 numpad encoding:
```
7(↖) 8(↑) 9(↗)
4(←) 5(·) 6(→)
1(↙) 2(↓) 3(↘)
```

Axis convention: dim1 = x (East+), dim2 = y (North+).
The table equals `OffsetArray(reshape(1:9, 3, 3), -2, -2)`:
```
# 1 4 7    W
# 2 5 8  S   N
# 3 6 9    E
```
"""
struct LDD <: FlowDirectionConvention end

_arrow(::Type{LDD}, d::Integer) = ('↖', '←', '↙', '↑', '·', '↓', '↗', '→', '↘')[d]

const _ldd_ci2dir = OffsetArray(reshape(UInt8.(1:9), 3, 3), -2, -2)
const _ldd_dir2ci = (
    CartesianIndex(-1, -1),  # 1 = SW
    CartesianIndex(0, -1),   # 2 = S
    CartesianIndex(1, -1),   # 3 = SE
    CartesianIndex(-1, 0),   # 4 = W
    CartesianIndex(0, 0),    # 5 = pit
    CartesianIndex(1, 0),    # 6 = E
    CartesianIndex(-1, 1),   # 7 = NW
    CartesianIndex(0, 1),    # 8 = N
    CartesianIndex(1, 1),    # 9 = NE
)

_ci_to_dir(::Type{LDD}, ci::CartesianIndex{2}) = _ldd_ci2dir[ci]
_dir_to_ci(::Type{LDD}, d::Integer) = _ldd_dir2ci[d]

"""
    D8D <: FlowDirectionConvention

D8 flow direction convention using power-of-2 encoding, clockwise from East:
```
 32(↖) 64(↑) 128(↗)
 16(←)  0(·)   1(→)
  8(↙)  4(↓)   2(↘)
```

Axis convention: dim1 = x (East+), dim2 = y (North+).
```
# 8 16  32    W
# 4  0  64  S   N
# 2  1 128    E
```

Values can be combined with bitwise OR to represent multiple flow directions
(e.g., `1 | 2 | 4 = 7` for E+SE+S). Use [`decompose`](@ref) to extract individual directions.
"""
struct D8D <: FlowDirectionConvention end

const _d8_single_arrows = ('↓', '↙', '←', '↖', '↑', '↗', '→', '↘')
const _d8_double_arrows = ('⇓', '⇙', '⇐', '⇖', '⇑', '⇗', '⇒', '⇘')

function _arrow(::Type{D8D}, d::Integer)
    d == 0 && return '·'
    ispow2(d) && return _d8_single_arrows[trailing_zeros(d) + 1]
    n = count_ones(d)
    if n == 2
        # Opposite pairs: double-headed arrows
        d == 1 | 16 && return '↔'  # E+W
        d == 4 | 64 && return '↕'  # S+N
        d == 2 | 32 && return '⤢'  # SE+NW
        d == 8 | 128 && return '⤡'  # SW+NE
        # Other pairs: double-stroke arrow for circular average direction
        b1 = trailing_zeros(d)
        b2 = trailing_zeros(d ⊻ (one(d) << b1))
        a1, a2 = b1 * 45.0, b2 * 45.0
        avg = atand(sind(a1) + sind(a2), cosd(a1) + cosd(a2))
        idx = mod(round(Int, mod(avg, 360) / 45, RoundNearestTiesUp), 8)
        return _d8_double_arrows[idx + 1]
    end
    return '✳'  # 3+ directions
end

const _d8_ci2dir = OffsetArray(UInt8[8 16 32; 4 0 64; 2 1 128], -2, -2)
const _d8_offsets = (
    CartesianIndex(1, 0),    # bit 0: 1 = E
    CartesianIndex(1, -1),   # bit 1: 2 = SE
    CartesianIndex(0, -1),   # bit 2: 4 = S
    CartesianIndex(-1, -1),  # bit 3: 8 = SW
    CartesianIndex(-1, 0),   # bit 4: 16 = W
    CartesianIndex(-1, 1),   # bit 5: 32 = NW
    CartesianIndex(0, 1),    # bit 6: 64 = N
    CartesianIndex(1, 1),    # bit 7: 128 = NE
)

_ci_to_dir(::Type{D8D}, ci::CartesianIndex{2}) = _d8_ci2dir[ci]
function _dir_to_ci(::Type{D8D}, d::Integer)
    d == 0 && return CartesianIndex(0, 0)
    ispow2(d) || throw(
        ArgumentError(
            "Combined direction $d cannot be converted to a single CartesianIndex, use decompose(D8D, $d)",
        ),
    )
    return _d8_offsets[trailing_zeros(d) + 1]
end

ismulti(::Type{D8D}) = true

"""Whether this direction represents a pit (no outflow)."""
ispit(d::FlowDirection{LDD}) = d.value == _ldd_ci2dir[CartesianIndex(0, 0)]
ispit(d::FlowDirection{D8D}) = iszero(d.value)

"""Whether this direction encodes exactly one flow direction (or pit)."""
issingle(d::FlowDirection{C}) where {C} = !ismulti(C) || ispow2(Int(d.value)) || ispit(d)

"""Number of flow directions encoded (0 for pit)."""
ndirections(d::FlowDirection{LDD}) = ispit(d) ? 0 : 1
ndirections(d::FlowDirection{D8D}) = count_ones(d.value)

"""
    decompose(d::Direction)
    decompose(R, d::Direction)
    decompose(R, d::Direction, center)

Lazily iterate over the individual single-direction values encoded by a direction,
or over relative grid indices of type `R`. The `center` argument supports grids
whose relative indices depend on the cell from which they originate.
For [`LDD`](@ref), returns a 1-tuple. For [`D8D`](@ref), iterates over each set bit.

# Examples
```julia
decompose(FlowDirection{D8D}(7))   # (→, ↘, ↓) → E + SE + S
decompose(FlowDirection{LDD}(9))   # (↘,) → SE
decompose(CartesianIndex{2}, FlowDirection{D8D}(3))
decompose(CartesianIndex{2}, FlowDirection{D8D}(3), CartesianIndex(2, 2))
```
"""
decompose(d::FlowDirection{LDD}) = (d,)
decompose(d::FlowDirection{D8D}) =
    (FlowDirection{D8D}(bit) for bit in decompose(D8D, d.value))
decompose(::Type{CartesianIndex{2}}, d::FlowDirection) =
    (CartesianIndex(direction) for direction in decompose(d))
decompose(::Type{R}, d::FlowDirection, center) where {R} = decompose(R, d)

struct DirectionBits{T <: Integer}
    value::T
end

Base.eltype(::Type{DirectionBits{T}}) where {T} = T
Base.length(directions::DirectionBits) = max(1, count_ones(directions.value))

function Base.iterate(directions::DirectionBits)
    iszero(directions.value) && return zero(directions.value), zero(directions.value)
    return _iterate_direction_bits(directions.value)
end
Base.iterate(::DirectionBits, remaining) =
    iszero(remaining) ? nothing : _iterate_direction_bits(remaining)

function _iterate_direction_bits(remaining)
    bit = one(remaining) << trailing_zeros(remaining)
    return bit, remaining ⊻ bit
end

decompose(::Type{D8D}, d::Integer) = DirectionBits(d)
