"""
    mask = skb(A; mean=mean(A))

Applies skewness balancing by [Bartels e.a (2006)](@cite bartelsDTMGenerationLIDAR2006) to `A`.
Improved the performance by applying a binary search to find the threshold value.

# Output
- `mask::BitMatrix` Mask of allowed values
"""
function skb(iA::AbstractArray; mean = _mean(iA))

    # Replace infinite values with maxintfloat
    mask = .!isfinite.(iA)
    if sum(mask) > 0
        A = copy(iA)
        A[mask] .= maxintfloat(eltype(A))
    else
        A = iA
    end

    # Sort A and get the indices
    I = sortperm(vec(A))
    AA = A[I]
    II = invperm(I)

    skew = 1
    splitby = 2
    len = step = i = length(AA) - sum(mask)

    # Search for the threshold value using a binary search
    while step >= 1
        skew = skewness(view(AA, firstindex(AA):i))
        step = len ÷ splitby
        splitby <<= 1
        if skew > 0
            i -= step
        else
            i += step
        end
        1 <= i <= len || break
    end
    if skew <= 0 || i < 1
        i += 1
    end

    fill!(mask, true)
    mask[I[i:end]] .= false
    return mask
end

function skb2(iA::AbstractArray; mean = mean(iA))
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
    fill!(m, true)
    m[I[i:end]] .= false
    return m
end

_mean(A::AbstractArray) = mean(filter(isfinite, vec(A)))

"""
    mask = skbr(A; iterations=10)

Applies recursive skewness balancing by [Bartels e.a (2010)](@cite bartelsThresholdfreeObjectGround2010) to `A`.
Applies `skb` `iterations` times to the object (non-terrain) mask, as to include
more (sloped) terrain.

# Output
- `mask::BitMatrix` Mask of allowed values
"""
function skbr(A::AbstractMatrix{<:Real}; iterations = 1, mean = _mean(A))
    @info mean
    terrain_mask = skb(A; mean)
    object_mask = .!terrain_mask
    while iterations > 1 && sum(object_mask) > 0
        @info "Iteration $iterations"
        terrain_mask[object_mask] .|= skb(A[object_mask])
        object_mask .= .!terrain_mask
        iterations -= 1
    end
    terrain_mask
end
