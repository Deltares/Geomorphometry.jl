const sqrt2 = sqrt(2.0)
const distance_8 = @SMatrix[sqrt2 1 sqrt2; 1 Inf 1; sqrt2 1 sqrt2]
const distance_4 = @SMatrix[Inf 1 Inf; 1 Inf 1; Inf 1 Inf]
const Δ = CartesianIndex(1, 1)

"Spread algorithms."
abstract type SpreadMethod end

Base.@kwdef struct Tomlin <: SpreadMethod end

"""
    spread(::Tomlin, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, method=Tomlin())

Total friction distance spread from `points` as by [Tomlin (1983)](@cite tomlin1983digital).
This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Matrix{<:Real}` Input Array
- `initial::Matrix{<:Real}` Initial values of the result
- `friction::Matrix{<:Real}` Resolution of cell size
- `res=1` Resolution or cell size
- `limit=Inf` Initial fill value
"""
function spread(
    ::Tomlin,
    locations::Vector{CartesianIndex{2}},
    initial::AbstractMatrix{T},
    friction::AbstractMatrix{<:Real};
    cellsize = cellsize(friction),
    limit = typemax(T),
) where {T <: Real}
    result = similar(friction)
    fill!(result, limit)
    result[locations] .= initial[locations]
    zone = zeros(Int32, size(friction))
    mask = falses(size(friction))

    II = CartesianIndices(size(friction))
    pq = FastPriorityQueue{eltype(friction)}(size(friction)...)
    for I in II[locations]
        enqueue!(pq, Tuple(I), result[I])
        mask[I] = true
    end

    δx, δy = abs.(cellsize)
    δxy = sqrt(δx^2 + δy^2)
    distances = @SMatrix[
        δxy δx δxy
        δy Inf δy
        δxy δx δxy
    ]

    while !isempty(pq)
        spread!(pq, mask, result, friction, zone, distances)
    end
    result
end

function spread!(pq, mask, result, friction, zone, distances)
    I = CartesianIndices(pq.index)[dequeue!(pq)]
    # I = dequeue!(pq)
    mask[I] = true
    patch = (I - Δ):(I + Δ)

    # New distance is cell_distance + average friction values
    for (li, i) in enumerate(patch)
        i in CartesianIndices(result) || continue
        mask[i] && continue
        nr = muladd(friction[i] + friction[I], distances[li] / 2, result[I])
        if nr < result[i]  # cells where new distance is lower
            result[i] = nr
            zone[i] = zone[I]
            haskey(pq, Tuple(i)) || enqueue!(pq, Tuple(i), result[i])
        end
    end
end

Base.@kwdef struct Eastman <: SpreadMethod
    iterations::Int = 3
end

"""
    spread(::Eastman, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, iterations=3)

Pushbroom method for friction costs as discussed by [Eastman (1989)](@cite eastman1989pushbroom).
This method should scale better (linearly) than the [Tomlin (1983)](@cite tomlin1983digital) method, but can require more
`iterations` than set by default (3) in the case of maze-like, uncrossable obstacles.

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Matrix{<:Real}` Input Array
- `initial::Matrix{<:Real}` Factor to exaggerate elevation
- `friction::Matrix{<:Real}` Resolution of cell size
- `res=1` Resolution or cell size
- `limit=Inf` Initial fill value
"""
function spread(
    e::Eastman,
    # points::AbstractMatrix{<:Real},
    locations::Vector{CartesianIndex{2}},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real};
    cellsize = cellsize(friction),
    limit = Inf,
)

    # ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ofriction[begin+1:end-1, begin+1:end-1] .= friction

    # result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    result = similar(friction)
    fill!(result, limit)

    # r = @view result[1:end-1, 1:end-1]
    # locations = points .> 0
    result[locations] .= initial[locations]

    # mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    # mask[begin+1:end-1, begin+1:end-1] .= false
    # mask = falses(size(points))

    # zone = OffsetMatrix(fill(0, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ozone = @view zone[1:end-1, 1:end-1]
    # ozone .= points
    # zone = copy(points)
    zone = zeros(Int32, size(friction))

    δx, δy = abs.(cellsize)
    δxy = sqrt(δx^2 + δy^2)
    distances = @SMatrix[
        δxy δx δxy
        δy Inf δy
        δxy δx δxy
    ]

    # minval, minidx = [0.0], [CartesianIndex(1, 1)]
    # x = @MMatrix zeros(3, 3)
    indices = CartesianIndices(size(friction))
    counts = 1
    for i in 1:(e.iterations)
        if iszero(counts)
            break
        end
        counts -= counts
        II = (i % 2 == 1) ? indices : reverse(indices)
        for I in II
            patch = (I - Δ):(I + Δ)
            for (li, i) in enumerate(patch)
                i in CartesianIndices(result) || continue
                nr = muladd(friction[i] + friction[I], distances[li] / 2, result[I])
                if nr < result[i]  # cells where new distance is lower
                    counts += 1
                    result[i] = nr
                    zone[i] = zone[I]
                end
            end
        end
    end
    result
end

"""
    spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Real; distance=Euclidean(), res=1.0)
    spread(points::Matrix{<:Real}, initial::Real, friction::Real; distance=Euclidean(), res=1.0)

Optimized (and more accurate) function based on the same friction everywhere.

When the friction is the same everywhere, there's no need for searching the shortest cost path,
as one can just take a direct line to the input points.

The calculated cost is more accurate, as there's no 'zigzag' from cell center to cell center.
"""
function spread(
    points::AbstractMatrix{<:Real},
    initial::AbstractMatrix{<:Real},
    friction::Real;
    distance = Euclidean(),
    cellsize = cellsize(friction),
    kwargs...,
)
    locations = findall(>(0), points)
    I = CartesianIndices(size(points))

    result = fill(Inf, size(points))
    for location in I[locations]
        for cell in I
            result[cell] = min(
                evaluate(distance, location.I, cell.I) * first(abs(cellsize)) * friction + initial[location],
                result[cell],
            )
        end
    end
    m = .~isfinite.(points)
    result[m] .= points[m]
    return result
end

Base.@kwdef struct FastSweeping <: SpreadMethod
    eps::Float64 = 1e-6
    debug::Bool = false
    iterations::Int = typemax(Int)
end

"""
    spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; cellsize=(1,1), limit=Inf, method=Tomlin())

Total friction distance spread from `points` from `initial` with `friction`.
By default uses Tomlin, see SpreadMethod for other algorithms.
"""
function spread(
    points::Vector{CartesianIndex{2}},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real};
    cellsize = cellsize(friction),
    limit = Inf,
    method = Tomlin(),
)
    spread(method, points, initial, friction; cellsize, limit)
end

function spread(
    points::AbstractMatrix{<:Real},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real};
    kwargs...,
)
    I = findall(>(0), points)
    return spread(I, initial, friction; kwargs...)
end

function spread(
    points::AbstractMatrix{<:Real},
    initial::Real,
    friction::AbstractMatrix{<:Real};
    kwargs...,
)
    init = Fill(initial, size(points))
    return spread(points, init, friction; kwargs...)
end

# function spread(
#     points::AbstractMatrix{<:Real},
#     initial::AbstractMatrix{<:Real},
#     friction::Real;
#     kwargs...,
# )
#     friction = Fill(friction, size(points))
#     return spread(points, initial, friction; kwargs...)
# end

function spread(points::AbstractMatrix{<:Real}, initial::Real, friction::Real; kwargs...)
    init = Fill(initial, size(points))
    fric = Fill(friction, size(points))
    return spread(points, init, fric; kwargs...)
end
