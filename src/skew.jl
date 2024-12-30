"""
    mask = skb(A; mean=mean(A))

Applies skewness balancing by *Bartels e.a (2006)* [^bartels2006] to `A`.
Improved the performance by applying a binary search to find the threshold value.

# Output
- `mask::BitMatrix` Mask of allowed values

[^bartels2006]: Bartels, M., Hong Wei, and D.C. Mason. 2006. “DTM Generation from LIDAR Data Using Skewness Balancing.” In 18th International Conference on Pattern Recognition (ICPR’06), 1:566–69. https://doi.org/10/cwk4v2.
"""
function skb(iA::AbstractArray; mean=mean(iA))
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
    len = length(AA) - sum(m)
    step = length(AA) - sum(m)
    i = length(AA) - sum(m)
    @info i

    while step >= 1
        s = skewness(AA[begin:i], mean)
        d <<= 1
        step = len ÷ d
        @info s, i, d
        if s > 0
            i -= step
        else
            i += step
        end
        1 <= i <= len || break
    end
    if s <= 0
        i += 1
    end

    fill!(m, true)
    m[I[i:end]] .= false
    return m
end

function skb2(iA::AbstractArray; mean=mean(iA))
    m = .!isfinite.(iA)
    if sum(m) > 0
        A = copy(iA)
        A[m] .= maxintfloat(eltype(A))
    else
        A = iA
    end
    I = sortperm(vec(A))
    AA = A[I]
    AAA = copy(AA)
    s = 1
    d = 2
    step = length(AA)

    S = sum(AA)
    N = length(AA)

    i = length(AA) - sum(m)
    while i > 0
        s = skewness(AA[begin:i])
        if s <= 0
            break
        end
        i -= 1
    end
    @info i
    fill!(m, true)
    m[I[i:end]] .= false
    return m
end


"""
    mask = skbr(A; iterations=10)

Applies recursive skewness balancing by *Bartels e.a (2006)* [^bartels2006] to `A`.
Applies `skb` `iterations` times to the object (non-terrain) mask, as to include
more (sloped) terrain.

# Output
- `mask::BitMatrix` Mask of allowed values

[^bartels2006]: Bartels, M., Hong Wei, and D.C. Mason. 2006. “DTM Generation from LIDAR Data Using Skewness Balancing.” In 18th International Conference on Pattern Recognition (ICPR’06), 1:566–69. https://doi.org/10/cwk4v2.
"""
function skbr(A; iterations=10, mean=mean(A))
    terrain_mask = skb(A; mean=mean)
    object_mask = .!terrain_mask
    while iterations > 1 && sum(object_mask) > 0
        # @info "$(round(Int, sum(object_mask) / length(object_mask) * 100))% objects..."
        terrain_mask[object_mask] = terrain_mask[object_mask] .| skb(A[object_mask])
        object_mask .= .!terrain_mask
        iterations -= 1
    end
    terrain_mask
end
