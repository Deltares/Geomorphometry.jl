const neib_8 = @SMatrix [1.0 1 1; 1 0 1; 1 1 1]
const neib_8_inf = @SMatrix [1.0 1 1; 1 Inf 1; 1 1 1]
const neib_8_mask = @SMatrix Bool[1 1 1; 1 0 1; 1 1 1]
const neib_8_dist = @SMatrix [√2 1 √2; 1 0 1; √2 1 √2]

const nbkernel = LocalFilters.Kernel{Int8, 2}(reshape(1:9, 3, 3))

abstract type DerivativeMethod end

"""Second order finite difference estimator using all 4 neighbors (Zevenbergen and Thorne, 1987)."""
struct ZevenbergenThorne <: DerivativeMethod end

"""Third order finite difference estimator using all 8 neighbors (Horn, 1981)."""
struct Horn <: DerivativeMethod end

"""Maximum Downward Gradient"""
struct MaximumDownwardGradient <: DerivativeMethod end
const MDG = MaximumDownwardGradient

struct LandSerf <: DerivativeMethod end

"""
    slope(dem::Matrix{<:Real}; cellsize=1.0, method=Horn(), exaggeration=1.0)

Slope is the rate of change between a cell and its neighbors as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function slope(
    dem::AbstractMatrix{<:Real};
    cellsize = 1.0,
    method::DerivativeMethod = Horn(),
    exaggeration = 1.0,
    xyratio=xyratio(dem)
)
    dst = copy(dem)
    slope!(method, dst, dem, cellsize, exaggeration, xyratio)
end

function slope!(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize, exaggeration, xyratio)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] = atand(
        √(((v[1] - v[2]) / (8 * v[5]))^2 + ((v[3] - v[4]) / (8 * v[5]))^2) * exaggeration,
    )
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function slope!(
    ::ZevenbergenThorne,
    dst,
    dem::AbstractMatrix{<:Real},
    cellsize,
    exaggeration,
    xyratio
)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    store!(d, i, v) = @inbounds d[i] =
        atand(√((v[1] / (2 * v[3]))^2 + (v[2] / (2 * v[3]))^2) * exaggeration)
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

function slope!(
    ::MaximumDownwardGradient,
    dst,
    dem::AbstractMatrix{<:Real},
    cellsize,
    exaggeration,
    xyratio
)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        m =
            max(
                abs(v[1]) / (v[5] * sqrt2),
                abs(v[2]) / v[5],
                abs(v[3]) / (v[5] * sqrt2),
                abs(v[4]) / v[5],
            ) / 2
        d[i] = atand(√(m^2 + m^2) * exaggeration)
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end

"""
    aspect(dem::Matrix{<:Real}, method=Horn())

Aspect is direction of [`slope`](@ref), as defined in Burrough, P. A., and McDonell, R. A., (1998, Principles of Geographical Information Systems).
"""
function aspect(dem::AbstractMatrix{<:Real}; method::DerivativeMethod = Horn(), xyratio=xyratio(dem))
    dst = copy(dem)
    aspect!(method, dst, dem, 1, xyratio)
end

# Useless, as there's no x/y component.
function aspect!(::MaximumDownwardGradient, dst, dem::AbstractMatrix{<:Real}, cellsize, xyratio)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx =
            max(
                abs(v[1]) / (v[5] * sqrt2),
                abs(v[2]) / v[5],
                abs(v[3]) / (v[5] * sqrt2),
                abs(v[4]) / v[5],
            ) / 2
        d[i] = compass(atand(δzδx, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, mdg, store!)
end

function aspect!(::Horn, dst, dem::AbstractMatrix{<:Real}, cellsize, xyratio)
    initial(A) =
        (zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = (v[1] - v[2]) / (8 * v[5])
        δzδy = (v[3] - v[4]) / (8 * v[5])
        d[i] = compass(atand(-δzδx, δzδy))
    end
    return localfilter!(dst, dem, nbkernel, initial, horn, store!)
end

function aspect!(::ZevenbergenThorne, dst, dem::AbstractMatrix{<:Real}, cellsize, xyratio)
    initial(A) = (zero(eltype(A)), zero(eltype(A)), cellsize)
    function store!(d, i, v)
        δzδx = v[1] / 2
        δzδy = v[2] / 2
        d[i] = compass(atand(δzδy, -δzδx))
    end
    return localfilter!(dst, dem, nbkernel, initial, zevenbergenthorne, store!)
end

@inline @inbounds function zevenbergenthorne(v, a, b)
    if b == 2
        return (v[1], v[2] + a, v[3])
    elseif b == 4
        return (v[1] + a, v[2], v[3])
    elseif b == 6
        return (v[1] - a, v[2], v[3])
    elseif b == 8
        return (v[1], v[2] - a, v[3])
    else
        return v
    end
end

@inline @inbounds function horn(v, a, b)
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
    if b == 1
        return (v[1] + a, v[2], v[3], v[4], v[5])
    elseif b == 2
        return (v[1], v[2], v[3], v[4] + a, v[5])
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
        return (v[1], v[2], v[3], v[4] - a, v[5])
    elseif b == 9
        return (v[1] - a, v[2], v[3], v[4], v[5])
    end
end

function compass(aspect)
    a = 90 - aspect
    a < 0 && (a += 360)
    return a
end

function aspect(compass::Real)
    return (450 - compass) % 360
end

@deprecate curvature(args...; kwargs...) laplacian(args...; gis=true, kwargs...)

function laplacian(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1, gis=false)
    mapstencil(x -> -2(_D(ZevenbergenThorne(),x,cellsize,cellsize) + _E(ZevenbergenThorne(),x,cellsize,cellsize)) * (gis ? 100 : 1), fitns(radius), dem)
end

fitns(R::Int) = NamedStencil(;
    Z1 = (-1R, -1R),
    Z2 = (-1R, 0),
    Z3 = (-1R, 1R),
    Z4 = (0, -1R),
    Z5 = (0, 0),
    Z6 = (0, 1R),
    Z7 = (1R, -1R),
    Z8 = (1R, 0),
    Z9 = (1R, 1R),
)


# https://www.spatialanalysisonline.com/HTML/index.html?profiles_and_curvature.htm
# https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-surface-parameters-works.htm

# Zevenbergen and Thorne (1987) method
# Surface fit by $Z = Ax²y² + Bx²y + Cxy² + Dx² + Ey² + Fxy + Gx + Hy + I$
_A(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) =
    ((s.Z1 + s.Z3 + s.Z7 + s.Z9) / 4 - (s.Z2 + s.Z4 + s.Z6 + s.Z8) / 2 + s.Z5) / δy^2*δx^2 
_B(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((s.Z1 + s.Z3 - s.Z7 - s.Z9) / 4 - (s.Z2 - s.Z8) / 2) / (δx^1.5 * δy^1.5)
_C(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((-s.Z1 + s.Z3 - s.Z7 + s.Z9) / 4 + (s.Z4 - s.Z6) / 2) / (δx ^1.5 * δy^1.5)
_D(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((s.Z4 + s.Z6) / 2 - s.Z5) / δx^2  # δ²z/δx²
_E(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = ((s.Z2 + s.Z8) / 2 - s.Z5) / δy^2  # δ²z/δy²
_F(::ZevenbergenThorne, s::Stencil, δx = 1, δy = 1) = (-s.Z1 + s.Z3 + s.Z7 - s.Z9) / (4δx * δy)   # δ²z/δxδy
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

coefficients(::LandSerf, s::Stencil, δx = 1, δy = 1) = (;
    a=_A(LandSerf(), s, δx, δy),
    b=_B(LandSerf(), s, δx, δy),
    c=_C(LandSerf(), s, δx, δy),
    d=_D(LandSerf(), s, δx, δy),
    e=_E(LandSerf(), s, δx, δy),
    f=_F(LandSerf(), s),
)

function profile(a, b, c, d, e)
    (-2 * (a * d^2 + c * d * e + b * e^2)) / ((d^2 + e^2) * (1 + d^2 + e^2)^1.5)
end

function _profile(s, cellsize)
    a, b, c, d, e = coefficients(LandSerf(), s, cellsize)
    profile(a, b, c, d, e)
end

"""
    profile_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)

Calculate normal slope line curvature (profile curvature)[^minar2020].

[^minar2020]: Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
"""
function profile_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)
    mapstencil(x->_profile(x, cellsize), fitns(radius), dem)
end

function tangential(a, b, c, d, e)
    -2 * (a * (e^2) - c * d * e + b * d^2) / ((d^2 + e^2) * sqrt(1 + d^2 + e^2))
end

function _tangential(s, cellsize)
    a, b, c, d, e = coefficients(LandSerf(), s, cellsize)
    tangential(a, b, c, d, e)
end

"""
    tangential_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)

Calculate normal contour curvature (tangential curvature)[^minar2020].

[^minar2020]: Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
"""
function tangential_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)
    mapstencil(x->_tangential(x, cellsize), fitns(radius), dem)
end


function plan(a, b, c, d, e)
    (2 * (b * d^2 - c * d * e + a * e^2)) / ((1 + d^2 + e^2)^1.5)
end

function _plan(s, cellsize)
    a, b, c, d, e = coefficients(LandSerf(), s, cellsize)
    plan(a, b, c, d, e)
end

"""
    plan_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)

Calculate projected contour curvature (plan curvature)[^minar2020].

[^minar2020]: Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
"""
function plan_curvature(dem::AbstractMatrix{<:Real}; cellsize = 1.0, radius=1)
    mapstencil(x->_plan(x, cellsize), fitns(radius), dem)
end

# """Calculate contour geodesic torsion or twisting curvature
# Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
# """
# function tgc(a, b, c, d, e)
#     2 * d * e * (a - b) - c * (d^2 - e^2) / ((d^2 + e^2) * (1 + d^2 + e^2))
# end

# """Calculate mean curvature
# Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
# """
# function kmean(a, b, c, d, e)
#     -(a * (1 + e^2) - c * d * e + b * (1 + d^2)) / (2sqrt((1 + d^2 + e^2)^3))
# end

# """Calculate unsphericity curvature
# Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
# """
# function ku(a, b, c, d, e)
#     sqrt(
#         ((a * (1 + e^2) - c * d * e + b * (1 + d^2)) / (sqrt((1 + d^2 + e^2)^3)))^2 -
#         ((4 * a * b - c^2) / (1 + d^2 + e^2)^2),
#     )
# end

# """Calculate min curvature
# Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
# """
# function kmin(a, b, c, d, e)
#     kmean(a, b, c, d, e) - ku(a, b, c, d, e)
# end

# """Calculate max curvature
# Minár, J., Evans, I.S., Jenčo, M., 2020. A comprehensive system of definitions of land surface (topographic) curvatures, with implications for their application in geoscience modelling and prediction. Earth-Science Reviews 211, 103414. https://doi.org/10.1016/j.earscirev.2020.103414
# """
# function kmax(a, b, c, d, e)
#     kmean(a, b, c, d, e) + ku(a, b, c, d, e)
# end