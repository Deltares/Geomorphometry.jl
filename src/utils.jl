using OffsetArrays
using ImageFiltering

"""Apply the opening operation to `A` with window size `ω`."""
function opening(A::Array{T,2}, ω::Integer) where T<:Real
    A = mapwindow(minimum, A, (ω,ω))  # erosion
    A = mapwindow(maximum, A, (ω,ω))  # dilation
    A
end

function circmask(n::Integer)
    kern = falses(-n:n, -n:n)
    for I in CartesianIndices(kern)
        i, j = Tuple(I)
        kern[I] = i^2 + j^2 <= n^2
    end
    return kern.parent
end

function opening_circ(A::Array{T,2}, ω::Integer) where T<:Real
    m = circmask(ω>>1)
    A = mapwindow(x->minimum(x[m]), A, (ω,ω))  # erosion
    A = mapwindow(x->maximum(x[m]), A, (ω,ω))  # dilation
    A
end
