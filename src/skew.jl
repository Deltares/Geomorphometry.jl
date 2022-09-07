"""
    mask = skb(A, mean)
    mask = skb(A)

Applies skewness balancing by *Bartels e.a (2006)* [^bartels2006] to `A`.
Improved the performance by applying a binary search to find the threshold value.

# Output
- `mask::BitMatrix` Mask of allowed values

[^bartels2006]: Bartels, M., Hong Wei, and D.C. Mason. 2006. “DTM Generation from LIDAR Data Using Skewness Balancing.” In 18th International Conference on Pattern Recognition (ICPR’06), 1:566–69. https://doi.org/10/cwk4v2.
"""
function skb(iA::AbstractArray{T}, mean::T) where {T<:Real}
    m = .!isfinite.(iA)
    if sum(m) > 0
        A = copy(iA)
        A[m] .= maxintfloat(eltype(A))
    else
        A = iA
    end
    I = sortperm(vec(A))
    AA = A[I]

    s = 1
    d = 2
    step = length(AA)
    i = length(AA)

    while step >= 1
        d <<= 1
        step = length(AA) ÷ d
        s = skewness(AA[begin:i], mean)
        if s > 0
            i -= step
        else
            i += step
        end
        1 <= i <= length(AA) || break
    end
    if s <= 0
        i += 1
    end

    fill!(m, true)
    m[I[i:end]] .= false
    return m
end

function skb(A::AbstractArray{T}) where {T<:Real}
    return skb(A, mean(A))
end
