using OffsetArrays
using ImageFiltering
using Distances
using PaddedViews

"""Apply the opening operation to `A` with window size `ω`."""
function opening(A::Array{T,2}, ω::Integer) where {T<:Real}
    A = mapwindow(minimum, A, (ω, ω))  # erosion
    A = mapwindow(maximum, A, (ω, ω))  # dilation
    A
end

"""Apply the opening operation to `A` with window size `ω`."""
function opening!(A::Array{T,2}, ω::Integer, out::Array{T,2}) where {T<:Real}
    mapwindow!(minimum, A, ω, out)  # erosion
    mapwindow!(maximum, out, ω, A)  # dilation
    A
end

"""Apply the opening operation to `A` with window size `ω`."""
function opening_circ!(A::Array{T,2}, ω::Integer, out::Array{T,2}) where {T<:Real}
    mapwindowcirc!(minimum_mask, A, ω, out, Inf)  # erosion
    mapwindowcirc!(maximum_mask, out, ω, A, -Inf)  # dilation
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
    Threads.@threads for I in R
        patch = max(I_first, I - Δ):min(I_last, I + Δ)
        out[I] = f(view(img, patch))
    end
    out
end

function mapwindowcirc!(f, img, window, out, fill = Inf)
    R = CartesianIndices(img)
    Δ = CartesianIndex(ntuple(x -> window ÷ 2, ndims(img)))

    w, h = size(img)
    A = PaddedView(fill, img, (-Δ[1]+1:w+Δ[1], -Δ[2]+1:h+Δ[2]))

    patch = (-Δ):(Δ)
    m = euclidean.(Tuple.(patch), Ref((0, 0))) .<= Δ[1]
    Threads.@threads for I in R
        patch = (I-Δ):(I+Δ)
        out[I] = f(view(A, patch), m)
    end
    out
end

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
