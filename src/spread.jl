const sqrt2 = sqrt(2.0)
const distance_8 = @SMatrix[sqrt2 1 sqrt2; 1 Inf 1; sqrt2 1 sqrt2]
const distance_4 = @SMatrix[Inf 1 Inf; 1 Inf 1; Inf 1 Inf]
const Δ = CartesianIndex(1, 1)

"Spread algorithms."
abstract type SpreadMethod end

"""
    Tomlin()

Friction-distance [`spread`](@ref) method based on a priority queue (Dijkstra-like search)
as described by [Tomlin (1983)](@cite tomlin1983digital). This is the default method.
"""
Base.@kwdef struct Tomlin <: SpreadMethod end

"""
    spread(::Tomlin, locations::Vector{CartesianIndex{2}}, initial::AbstractMatrix{<:Real}, friction::AbstractMatrix{<:Real}; cellsize=cellsize(friction), limit=typemax(T))

Total friction distance spread from `locations` as described by [Tomlin (1983)](@cite tomlin1983digital).
This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).

# Output
- `Matrix{<:Real}` Total friction distance

# Arguments
- `locations::Vector{CartesianIndex{2}}` Source cells to spread from
- `initial::AbstractMatrix{<:Real}` Initial values of the result at the source cells
- `friction::AbstractMatrix{<:Real}` Friction (cost) per cell
- `cellsize=cellsize(friction)` Cell size, used to scale distances
- `limit=typemax(T)` Initial fill value for unreached cells
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

"""
    Eastman(; iterations=3)

Friction-distance [`spread`](@ref) method using the pushbroom approach of
[Eastman (1989)](@cite eastman1989pushbroom). Scales better (linearly) than [`Tomlin`](@ref),
but may require more `iterations` for maze-like, uncrossable obstacles.
"""
Base.@kwdef struct Eastman <: SpreadMethod
    iterations::Int = 3
end

"""
    spread(::Eastman, locations::Vector{CartesianIndex{2}}, initial::AbstractMatrix{<:Real}, friction::AbstractMatrix{<:Real}; cellsize=cellsize(friction), limit=Inf)

Pushbroom method for friction costs as discussed by [Eastman (1989)](@cite eastman1989pushbroom).
This method should scale better (linearly) than the [Tomlin (1983)](@cite tomlin1983digital) method, but can require more
`iterations` than set by default (3) in the case of maze-like, uncrossable obstacles.

# Output
- `Matrix{<:Real}` Total friction distance

# Arguments
- `locations::Vector{CartesianIndex{2}}` Source cells to spread from
- `initial::AbstractMatrix{<:Real}` Initial values of the result at the source cells
- `friction::AbstractMatrix{<:Real}` Friction (cost) per cell
- `cellsize=cellsize(friction)` Cell size, used to scale distances
- `limit=Inf` Initial fill value for unreached cells
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
    spread(points::AbstractMatrix{<:Real}, initial::AbstractMatrix{<:Real}, friction::Real; distance=Euclidean(), cellsize=cellsize(friction))
    spread(points::AbstractMatrix{<:Real}, initial::Real, friction::Real; distance=Euclidean(), cellsize=cellsize(friction))

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

"""
    FastSweeping(; eps=1e-6, debug=false, iterations=typemax(Int))

Friction-distance [`spread`](@ref) method using an iterative fast sweeping scheme. Sweeps
the grid in alternating directions until the result converges within `eps` or `iterations`
is reached.
"""
Base.@kwdef struct FastSweeping <: SpreadMethod
    eps::Float64 = 1e-6
    debug::Bool = false
    iterations::Int = typemax(Int)
end

"""
    spread(points, initial, friction::AbstractMatrix{<:Real}; cellsize=cellsize(friction), limit=Inf, method=Tomlin())

Total friction distance spread from `points` from `initial` with `friction`.

`points` may be a `Vector{CartesianIndex{2}}` of source cells, or an `AbstractMatrix{<:Real}`
in which case the cells with values greater than zero are taken as sources. `initial` is
either a matrix of source values or a single scalar applied to all sources. By default uses
[`Tomlin`](@ref); see `SpreadMethod` for other algorithms.
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
