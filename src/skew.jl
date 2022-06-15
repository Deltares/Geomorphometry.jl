# Skewness balancing as by Bartels (2006)
"""
```
mask = skb(A, mean)
mask = skb(A)
```
Applies skewness balancing by *Bartels e.a (2006)* [^bartels2006] to `A`.

# Output
- `mask::BitMatrix` Maximum allowable values

Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.

[^bartels2006]: Bartels, M., Hong Wei, and D.C. Mason. 2006. “DTM Generation from LIDAR Data Using Skewness Balancing.” In 18th International Conference on Pattern Recognition (ICPR’06), 1:566–69. https://doi.org/10/cwk4v2.
"""
function skb(A::AbstractArray{T}, mean::T) where {T<:Real}
    I = sortperm(vec(A))
    AA = A[I]
    s = 1
    for i ∈ length(A):-1:1
        X = skewness(AA[begin:i], mean)
        if X <= 0
            s = i + 1
            break
        end
    end
    m = trues(size(A))
    m[I[s+1:end]] .= false
    return m
end

function skb(A::AbstractArray{T}) where {T<:Real}
    return skb(A, mean(A))
end
