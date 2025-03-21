# const neib_8 = @SMatrix [1.0 1 1; 1 0 1; 1 1 1]
# const neib_8_inf = @SMatrix [1.0 1 1; 1 Inf 1; 1 1 1]
# const neib_8_mask = @SMatrix Bool[1 1 1; 1 0 1; 1 1 1]
const neib_8_dist = @SMatrix [√2 1 √2; 1 0 1; √2 1 √2]

const nbkernel = LocalFilters.Kernel{Int8, 2}(reshape(1:9, 3, 3))

abstract type DerivativeMethod end

"""Second order finite difference estimator using all 4 neighbors by [Zevenbergen and Thorne, (1987)](@cite zevenbergen1987quantitative)."""
struct ZevenbergenThorne <: DerivativeMethod end

"""Third order finite difference estimator using all 8 neighbors by [Horn, (1981)](@cite hornHillShadingReflectance1981)."""
struct Horn <: DerivativeMethod end

"""Maximum Downward Gradient"""
struct MaximumDownwardGradient <: DerivativeMethod end
const MDG = MaximumDownwardGradient

struct LandSerf <: DerivativeMethod end

struct GDAL <: DerivativeMethod end

"""
    slope(dem::Matrix{<:Real}; cellsize=cellsize(dem), method=Horn(), exaggeration=1.0)

Slope is the rate of change between a cell and its neighbors.
"""
function slope(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    method::DerivativeMethod = Horn(),
    exaggeration = 1.0,
    direction::Union{Nothing, Real} = nothing,
)
    dst = similar(dem, Float32)
    slope!(method, dst, dem, cellsize, exaggeration, direction)
end

function slope!(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize, exaggeration, direction)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)

    store!(d, i, v) = @inbounds d[i] =
        isnothing(direction) ?
        atand(
            √(((v[3] - v[4]) / (8 * v[5][1]))^2 + ((v[1] - v[2]) / (8 * v[5][2]))^2) *
            exaggeration,
        ) :
        atand(
            (
                -((v[3] - v[4]) / (8 * v[5][1])) * cosd(_aspect(direction)) +
                ((v[1] - v[2]) / (8 * v[5][2])) * sind(_aspect(direction))
            ) * exaggeration,
        )

    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function slope!(
    ::ZevenbergenThorne,
    dst,
    dem::AbstractMatrix{<:Real},
    cellsize,
    exaggeration,
    direction,
)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] =
        isnothing(direction) ?
        atand(√((v[1] / (2 * v[3][1]))^2 + (v[2] / (2 * v[3][2]))^2) * exaggeration) :
        atand(
            -(v[1] / (2 * v[3][1])) * cosd(_aspect(direction)) +
            (v[2] / (2 * v[3][2])) * sind(_aspect(direction)),
        ) * exaggeration
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

function slope!(
    ::MaximumDownwardGradient,
    dst,
    dem::AbstractMatrix{<:Real},
    cellsize,
    exaggeration,
    direction,
)
    isnothing(direction) || throw(ArgumentError("Direction is not supported for MDG."))
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        m =
            max(
                abs(v[1]) / sqrt(v[5][1]^2 * v[5][2]^2),  # \
                abs(v[2] / v[5][2]),  # |
                abs(v[3]) / sqrt(v[5][1]^2 * v[5][2]^2),  # /
                abs(v[4] / v[5][1]),  # -
            ) / 2
        d[i] = atand(√(m^2 + m^2) * exaggeration)
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end

"""
    aspect(dem::Matrix{<:Real}, method=Horn())

Aspect is direction of [`slope`](@ref).
"""
function aspect(
    dem::AbstractMatrix{<:Real};
    method::DerivativeMethod = Horn(),
    cellsize = cellsize(dem),
)
    dst = similar(dem, Float32)
    aspect!(method, dst, dem, cellsize)
end

# Useless, as there's no x/y component.
function aspect!(::MaximumDownwardGradient, dst, dem::AbstractMatrix{<:Real}, cellsize)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = -Inf
        aspect = NaN

        nbs = (
            (v[1] / sqrt(v[5][1]^2 * v[5][2]^2), 45), # /
            (v[2] / v[5][2], 90), # -
            (v[3] / sqrt(v[5][1]^2 * v[5][2]^2), 45 + 90), # \
            (v[4] / v[5][1], 180), # |
        )

        for (candidate, orientation) in nbs
            if abs(candidate) > δzδx
                δzδx = abs(candidate)
                aspect = candidate > 0 ? orientation : (orientation + 180) % 360
            end
        end
        d[i] = aspect
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end

function aspect!(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize)
    # 1north 2south 3east 4west 5cellsize
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = (v[3] - v[4]) / (8 * v[5][1])
        δzδy = (v[1] - v[2]) / (8 * v[5][2])
        d[i] = compass(atand(-δzδy, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function aspect!(::ZevenbergenThorne, dst, dem::AbstractMatrix{<:Real}, cellsize)
    # 1ew 2ns 3cellsize
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = v[1] / (2 * v[3][1])
        δzδy = v[2] / (2 * v[3][2])
        d[i] = compass(atand(-δzδy, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

@inline @inbounds function zevenbergenthorne(v, a, b)
    # 1 4 7    W
    # 2 5 8  S   N
    # 3 6 9    E
    if b == 2  # South
        return (v[1], v[2] - a, v[3])
    elseif b == 4  # West
        return (v[1] - a, v[2], v[3])
    elseif b == 6  # East
        return (v[1] + a, v[2], v[3])
    elseif b == 8  # North
        return (v[1], v[2] + a, v[3])
    else
        return v
    end
end

@inline @inbounds function horn(v, a, b)
    # 1 4 7    W
    # 2 5 8  S   N
    # 3 6 9    E
    if b == 1
        return (v[1], v[2] + a, v[3], v[4] + a, v[5])
    elseif b == 2
        return (v[1], v[2] + 2a, v[3], v[4], v[5])
    elseif b == 3
        return (v[1], v[2] + a, v[3] + a, v[4], v[5])
    elseif b == 4
        return (v[1], v[2], v[3], v[4] + 2a, v[5])
    elseif b == 5
        return v
    elseif b == 6
        return (v[1], v[2], v[3] + 2a, v[4], v[5])
    elseif b == 7
        return (v[1] + a, v[2], v[3], v[4] + a, v[5])
    elseif b == 8
        return (v[1] + 2a, v[2], v[3], v[4], v[5])
    elseif b == 9
        return (v[1] + a, v[2], v[3] + a, v[4], v[5])
    end
end

@inline @inbounds function mdg(v, a, b)
    # 1 4 7    W
    # 2 5 8  S   N
    # 3 6 9    E
    if b == 1
        return (v[1] + a, v[2], v[3], v[4], v[5])
    elseif b == 2
        return (v[1], v[2], v[3], v[4] - a, v[5])
    elseif b == 3
        return (v[1], v[2], v[3] - a, v[4], v[5])
    elseif b == 4
        return (v[1], v[2] + a, v[3], v[4], v[5])
    elseif b == 5
        return v
    elseif b == 6
        return (v[1], v[2] - a, v[3], v[4], v[5])
    elseif b == 7
        return (v[1], v[2], v[3] + a, v[4], v[5])
    elseif b == 8
        return (v[1], v[2], v[3], v[4] + a, v[5])
    elseif b == 9
        return (v[1] - a, v[2], v[3], v[4], v[5])
    end
end

function compass(aspect)
    a = 270 - aspect
    (a + 360) % 360
end

function _aspect(compass::Real)
    return (compass - 90) % 360
end

@deprecate curvature(args...; kwargs...) laplacian(args...; gis = true, kwargs...)

function laplacian(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    radius = 1,
    gis = false,
    direction = nothing,
)
    mapstencil(
        x ->
            -2(
                _D(ZevenbergenThorne(), x, cellsize[1], cellsize[2]) *
                (isnothing(direction) ? 1 : cos(_aspect(direction))) +
                _E(ZevenbergenThorne(), x, cellsize[1], cellsize[2]) *
                (isnothing(direction) ? 1 : sin(_aspect(direction)))
            ) * (gis ? 100 : 1),
        scaled8nb(radius),
        dem,
    )
end

scaled8nb(R::Int) = NamedStencil(;
    Z1 = (-1R, -1R),
    Z2 = (0, -1R),  # South
    Z3 = (1R, -1R),
    Z4 = (-1R, 0),  # West
    Z5 = (0, 0),  # Center
    Z6 = (1R, 0),  # East
    Z7 = (-1R, 1R),
    Z8 = (0, 1R),  # North
    Z9 = (1R, 1R),
)

# https://www.spatialanalysisonline.com/HTML/index.html?profiles_and_curvature.htm
# https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-surface-parameters-works.htm

# Zevenbergen and Thorne (1987) method
# Surface fit by $Z = Ax²y² + Bx²y + Cxy² + Dx² + Ey² + Fxy + Gx + Hy + I$
_A(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) =
    ((s.Z1 + s.Z3 + s.Z7 + s.Z9) / 4 - (s.Z2 + s.Z4 + s.Z6 + s.Z8) / 2 + s.Z5) / δy^2 * δx^2
_B(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) =
    ((s.Z1 + s.Z3 - s.Z7 - s.Z9) / 4 - (s.Z2 - s.Z8) / 2) / (δx^1.5 * δy^1.5)
_C(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) =
    ((-s.Z1 + s.Z3 - s.Z7 + s.Z9) / 4 + (s.Z4 - s.Z6) / 2) / (δx^1.5 * δy^1.5)
_D(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((s.Z4 + s.Z6) / 2 - s.Z5) / δx^2  # δ²z/δx²
_E(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((s.Z2 + s.Z8) / 2 - s.Z5) / δy^2  # δ²z/δy²
_F(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) =
    (-s.Z1 + s.Z3 + s.Z7 - s.Z9) / (4δx * δy)   # δ²z/δxδy
_G(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = (-s.Z4 + s.Z6) / 2δx  # δz/δx
_H(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = (s.Z2 - s.Z8) / 2δy  # δz/δy
_I(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = s.Z5

# Landserf 
# Surface fit by $ z = Ax² + By² + Cxy + Dx + Ey + F $
# Note that the approximations ABCDEF here are (almost) identical to the DEFGHI in Zevenbergen and Thorne
_A(::LandSerf, s::Stencil, δx = 1, δy = 1) = (s.Z6 + s.Z4 - 2s.Z5) / δx^2 # δ²z/δx²
_B(::LandSerf, s::Stencil, δx = 1, δy = 1) = (s.Z8 + s.Z2 - 2s.Z5) / δy^2 # δ²z/δy²
_C(::LandSerf, s::Stencil, δx = 1, δy = 1) = (s.Z9 - s.Z7 - s.Z3 + s.Z1) / (4δx * δy)  # δ²z/δxδy
_D(::LandSerf, s::Stencil, δx = 1, δy = 1) = (s.Z6 - s.Z4) / 2δx  # δz/δx
_E(::LandSerf, s::Stencil, δx = 1, δy = 1) = (s.Z8 - s.Z2) / 2δy  # δz/δy
_F(::LandSerf, s::Stencil) = s.Z5

function coefficients(::LandSerf, s::Stencil, δx = 1, δy = 1, direction = nothing)
    if isnothing(direction)
        (;
            a = _A(LandSerf(), s, δx, δy),
            b = _B(LandSerf(), s, δx, δy),
            c = _C(LandSerf(), s, δx, δy),
            d = _D(LandSerf(), s, δx, δy),
            e = _E(LandSerf(), s, δx, δy),
            f = _F(LandSerf(), s),
        )
    else
        direction = _aspect(direction)
        (;
            a = _A(LandSerf(), s, δx, δy) * cos(direction)^2,
            b = _B(LandSerf(), s, δx, δy) * sin(direction)^2,
            c = _C(LandSerf(), s, δx, δy) * sin(direction) * cos(direction),
            d = _D(LandSerf(), s, δx, δy) * cos(direction),
            e = _E(LandSerf(), s, δx, δy) * sin(direction),
            f = _F(LandSerf(), s),
        )
    end
end

function profile(a, b, c, d, e)
    (-2 * (a * d^2 + c * d * e + b * e^2)) / ((d^2 + e^2) * (1 + d^2 + e^2)^1.5)
end

function _profile(s, δx, δy, direction)
    a, b, c, d, e = coefficients(LandSerf(), s, δx, δy, direction)
    profile(a, b, c, d, e)
end

"""
    profile_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)

Calculate normal slope line curvature (profile curvature) as defined by [Minár et al., (2020)](@cite minarComprehensiveSystemDefinitions2020).
"""
function profile_curvature(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    radius = 1,
    direction = nothing,
)
    mapstencil(
        x -> _profile(x, cellsize[1], cellsize[2], direction),
        scaled8nb(radius),
        dem,
    )
end

function tangential(a, b, c, d, e)
    -2 * (a * (e^2) - c * d * e + b * d^2) / ((d^2 + e^2) * sqrt(1 + d^2 + e^2))
end

function _tangential(s, δx, δy, direction)
    a, b, c, d, e = coefficients(LandSerf(), s, δx, δy, direction)
    tangential(a, b, c, d, e)
end

"""
    tangential_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)

Calculate normal contour curvature (tangential curvature) as defined by [Minár et al., (2020)](@cite minarComprehensiveSystemDefinitions2020).
"""
function tangential_curvature(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    radius = 1,
    direction = nothing,
)
    mapstencil(
        x -> _tangential(x, cellsize[1], cellsize[2], direction),
        scaled8nb(radius),
        dem,
    )
end

function plan(a, b, c, d, e)
    (2 * (b * d^2 - c * d * e + a * e^2)) / ((1 + d^2 + e^2)^1.5)
end

function _plan(s, δx, δy, direction)
    a, b, c, d, e = coefficients(LandSerf(), s, δx, δy, direction)
    plan(a, b, c, d, e)
end

"""
    plan_curvature(dem::AbstractMatrix{<:Real}; cellsize = cellsize(dem), radius=1)

Calculate projected contour curvature (plan curvature) as defined by [Minár et al., (2020)](@cite minarComprehensiveSystemDefinitions2020).
"""
function plan_curvature(
    dem::AbstractMatrix{<:Real};
    cellsize = cellsize(dem),
    radius = 1,
    direction = nothing,
)
    mapstencil(x -> _plan(x, cellsize[1], cellsize[2], direction), scaled8nb(radius), dem)
end
