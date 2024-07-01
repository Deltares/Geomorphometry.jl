"""Apply the opening operation to `A` with window size `ω`."""
function opening(A::AbstractArray{T,2}, ω::Integer) where {T<:Real}
    A = mapwindow(minimum, A, (ω, ω))  # erosion
    A = mapwindow(maximum, A, (ω, ω))  # dilation
    A
end

"""Apply the opening operation to `A` with window size `ω`."""
function opening!(A::AbstractArray{T,2}, ω::Integer, out::AbstractArray{T,2}) where {T<:Real}
    mapwindow_sep!(minimum, A, ω, out, Inf)  # erosion
    mapwindow_sep!(maximum, out, ω, A, -Inf)  # dilation
    A
end

"""Apply the opening operation to `A` with window size `ω`."""
function opening_circ!(A::AbstractArray{T,2}, ω::Integer, out::AbstractArray{T,2}) where {T<:Real}
    mapwindowcirc_approx!(minimum_mask, A, ω, out, Inf)  # erosion
    mapwindowcirc_approx!(maximum_mask, out, ω, A, -Inf)  # dilation
    A
end

function circmask(n::Integer)
    # TODO This could be precomputed for first N integers
    kern = falses(-n:n, -n:n)
    for I in CartesianIndices(kern)
        i, j = Tuple(I)
        kern[I] = i^2 + j^2 <= n^2
    end
    return kern.parent
end

function opening_circ(A::Array{T,2}, ω::Integer) where {T<:Real}
    m = circmask(ω >> 1)
    A = mapwindow(x -> minimum(x[m]), A, (ω, ω))  # erosion
    A = mapwindow(x -> maximum(x[m]), A, (ω, ω))  # dilation
    A
end


# First discussed here https://github.com/JuliaImages/ImageFiltering.jl/issues/179
function mapwindow!(f, img, window, out)
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(ntuple(x -> window ÷ 2, ndims(img)))
    @inbounds for I in R
        patch = max(I_first, I - Δ):min(I_last, I + Δ)
        out[I] = f(view(img, patch))
    end
    out
end

function mapwindow_stack!(f, img, window, out)
    R = CartesianIndices(img)
    I_first, I_last = first(R), last(R)
    Δ = CartesianIndex(1, 1)
    out2 = copy(img)
    iterations = window:-2:3
    @inbounds for _ in iterations  # repeat 3x3 window
        for I in R
            patch = max(I_first, I - Δ):min(I_last, I + Δ)
            out[I] = f(view(out2, patch))
        end
        out2 .= out
    end
    out
end

function mapwindow_sep!(f, img, window, out, fill=Inf)
    Δ = window ÷ 2

    w, h = size(img)
    A = PaddedView(fill, img, (-Δ+1:w+Δ, -Δ+1:h+Δ))
    out2 = copy(out)

    # Maximum/minimum is seperable into 1d
    @inbounds for i in 1:h, j in 1:w
        out2[j, i] = f(@view A[j-Δ:j+Δ, i])
    end
    A = PaddedView(fill, out2, (-Δ+1:w+Δ, -Δ+1:h+Δ))
    @inbounds for j in 1:w, i in 1:h
        out[j, i] = f(@view A[j, i-Δ:i+Δ])
    end
    out
end


function mapwindowcirc!(f, img, window, out, fill=Inf)
    R = CartesianIndices(img)
    Δ = CartesianIndex(ntuple(x -> window ÷ 2, ndims(img)))

    w, h = size(img)
    A = PaddedView(fill, img, (-Δ[1]+1:w+Δ[1], -Δ[2]+1:h+Δ[2]))

    m = euclidean.(Tuple.(-Δ:Δ), Ref((0, 0))) .<= Δ[1]
    @inbounds for I in R
        patch = (I-Δ):(I+Δ)
        out[I] = f(view(A, patch), m)
    end
    out
end


function mapwindowcirc_approx!(f, img, window, out, fill=Inf)
    R = CartesianIndices(img)
    Δ = CartesianIndex(1, 1)

    w, h = size(img)

    iterations = window:-2:3
    A = PaddedView(fill, img, (-Δ[1]+1:w+Δ[1], -Δ[2]+1:h+Δ[2]))

    m = euclidean.(Tuple.(-Δ:Δ), Ref((0, 0))) .<= Δ[1]
    @inbounds for _ in iterations  # repeat 3x3 window
        for I in R
            patch = (I-Δ):(I+Δ)
            out[I] = f(view(A, patch), m)
        end
        img .= out
    end
    out
end

# Functions for future changes, based on LocalFiltering
# function opening_circ_approx2!(A::Array{T,2}, ω::Integer, out::Array{T,2}) where {T<:Real}
#     iterations = ω:-2:3

#     B = circmask(1)
#     for _ in iterations
#         localfilter!(out, A, B,
#             (a) -> typemax(a),
#             (v, a, b) -> b ? min(v, a) : v,
#             (d, i, v) -> @inbounds(d[i] = v))
#         A, out = out, A
#     end
#     for _ in iterations
#         localfilter!(out, A, B,
#             (a) -> typemin(a),
#             (v, a, b) -> b ? max(v, a) : v,
#             (d, i, v) -> @inbounds(d[i] = v))
#         A, out = out, A
#     end
#     A
# end

# function opening_circ2!(A::Array{T,2}, ω::Integer, out::Array{T,2}) where {T<:Real}
#     B = circmask(ω >> 1)
#     localfilter!(out, A, B,
#         (a) -> typemax(a),
#         (v, a, b) -> b ? min(v, a) : v,
#         (d, i, v) -> @inbounds(d[i] = v))
#     localfilter!(A, out, B,
#         (a) -> typemin(a),
#         (v, a, b) -> b ? max(v, a) : v,
#         (d, i, v) -> @inbounds(d[i] = v))
#     A
# end

@inline function maximum_mask(x, m)
    o = -Inf
    @inbounds for I in eachindex(x)
        o = ifelse(m[I], max(o, x[I]), o)
    end
    o
end

@inline function minimum_mask(x, m)
    o = Inf
    @inbounds for I in eachindex(x)
        o = ifelse(m[I], min(o, x[I]), o)
    end
    o
end
