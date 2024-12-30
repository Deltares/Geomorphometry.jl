const sqrt2 = sqrt(2.0)
const distance_8 = @SMatrix[sqrt2 1 sqrt2; 1 Inf 1; sqrt2 1 sqrt2]
const distance_4 = @SMatrix[Inf 1 Inf; 1 Inf 1; Inf 1 Inf]
const Δ = CartesianIndex(1, 1)

"Spread algorithms."
abstract type SpreadMethod end

Base.@kwdef struct Tomlin <: SpreadMethod end

"""
    spread(::Tomlin, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, method=Tomlin())

Total friction distance spread from `points` as by *Tomlin (1983)* [^tomlin1983].
This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Matrix{<:Real}` Input Array
- `initial::Matrix{<:Real}` Initial values of the result
- `friction::Matrix{<:Real}` Resolution of cell size
- `res=1` Resolution or cell size
- `limit=Inf` Initial fill value

[^tomlin1983]: Tomlin, Charles Dana. 1983. Digital Cartographic Modeling Techniques in Environmental Planning. Yale University.
"""
function spread(
    t::Tomlin,
    locations::Vector{CartesianIndex{2}},
    initial::AbstractMatrix{T},
    friction::AbstractMatrix{<:Real};
    res = 1,
    limit = typemax(T),
) where {T <: Real}

    # ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ofriction[begin+1:end-1, begin+1:end-1] .= friction

    # result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    result = fill(limit, size(friction))
    # r = @view result[1:end-1, 1:end-1]
    # locations = points .> 0
    result[locations] .= initial[locations]

    # zone = OffsetMatrix(fill(0, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ozone = @view zone[1:end-1, 1:end-1]
    # ozone .= points
    zone = zeros(Int32, size(friction))

    # Construct stack for locations
    # mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    # mask[begin+1:end-1, begin+1:end-1] .= false
    mask = falses(size(friction))

    II = CartesianIndices(size(friction))
    # pq = PriorityQueue{CartesianIndex,eltype(friction)}()
    pq = FastPriorityQueue{eltype(friction)}(size(locations)...)
    # pq = PriorityQueue{CartesianIndex{2},eltype(friction)}()
    for I in II[locations]
        enqueue!(pq, Tuple(I), result[I])
        mask[I] = true
    end
    # Step 1: Set the distance of the starting node to 0 and the distances of all other nodes to the highest value possible.

    # Step 3: For each of the active node’s adjacent neighbors, set its distance to whichever is
    # less: its current distance value or the sum of the distance of the active node plus the
    # weight of the arc from the active node to that neighbor.
    while !isempty(pq)
        spread!(pq, mask, result, friction, zone, res)
    end
    result, zone
end

function spread!(pq, mask, result, friction, zone, res)
    I = CartesianIndices(pq.index)[dequeue!(pq)]
    # I = dequeue!(pq)
    mask[I] = true
    patch = (I - Δ):(I + Δ)

    # New distance is cell_distance + average friction values
    for (li, i) in enumerate(patch)
        i in CartesianIndices(result) || continue
        mask[i] && continue
        nr = muladd(friction[i] + friction[I], res / 2 * distance_8[li], result[I])
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

Pushbroom method for friction costs as discussed by *Eastman (1989) [^eastman1989].
This method should scale much better (linearly) than the [^tomlin1983] method, but can require more
`iterations` than set by default (3) in the case of maze-like, uncrossable obstacles.

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Matrix{<:Real}` Input Array
- `initial::Matrix{<:Real}` Factor to exaggerate elevation
- `friction::Matrix{<:Real}` Resolution of cell size
- `res=1` Resolution or cell size
- `limit=Inf` Initial fill value

[^eastman1989]: Eastman, J. Ronald. 1989. ‘Pushbroom Algorithms for Calculating Distances in Raster Grids’. In Proceedings, Autocarto, 9:288–97.
"""
function spread(
    ::Eastman,
    points::AbstractMatrix{<:Real},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real};
    res = 1,
    limit = Inf,
)

    # ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ofriction[begin+1:end-1, begin+1:end-1] .= friction

    # result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    result = fill(limit, size(friction))

    # r = @view result[1:end-1, 1:end-1]
    locations = points .> 0
    result[locations] .= initial[locations]

    # mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    # mask[begin+1:end-1, begin+1:end-1] .= false
    # mask = falses(size(points))

    # zone = OffsetMatrix(fill(0, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    # ozone = @view zone[1:end-1, 1:end-1]
    # ozone .= points
    zone = copy(points)
    # minval, minidx = [0.0], [CartesianIndex(1, 1)]
    # x = @MMatrix zeros(3, 3)
    indices = CartesianIndices(size(points))
    counts = 1
    for i in 1:iterations
        if iszero(counts)
            break
        end
        counts -= counts
        II = (i % 2 == 1) ? indices : reverse(indices)
        for I in II
            patch = (I - Δ):(I + Δ)
            for (li, i) in enumerate(patch)
                i in CartesianIndices(result) || continue
                nr = muladd(friction[i] + friction[I], res / 2 * distance_8[li], result[I])
                if nr < result[i]  # cells where new distance is lower
                    counts += 1
                    result[i] = nr
                    zone[i] = zone[I]
                end
            end
        end
    end
    result, zone
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
    res = 1.0,
    kwargs...,
)
    locations = findall(>(0), points)
    I = CartesianIndices(size(points))

    result = fill(Inf, size(points))
    for location in I[locations]
        for cell in I
            result[cell] = min(
                evaluate(distance, location.I, cell.I) * res * friction + initial[location],
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
    spread(::FastSweeping, points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf)

Total friction distance spread from `points` as by algorithm *Tomlin (1983)* [^tomlin1983].
This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).
See SpreadMethod for other algorithms.

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Vector{CartesianIndex}` Input Array
- `initial::AbstractVector{<:Real}` Initial values of the result
- `friction::Matrix{<:Real}` Friction map

[^hongkai20025]: Zhao, Hongkai (2005). "A fast sweeping method for Eikonal equations". Mathematics of Computation. 74 (250): 603–627. DOI: 10.1090/S0025-5718-04-01678-3
"""
function spread(fs::FastSweeping, points, initial, friction)
    solver = Eikonal.FastSweeping(friction)
    for I in points
        solver.t[I] = initial[I]
    end
    sweep!(solver; nsweeps = fs.iterations, fs.debug, fs.eps)
    @view solver.t[2:end, 2:end]
end

"""
    spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, method=Tomlin())

Total friction distance spread from `points` from `initial` with `friction`.
By default uses Tomlin, see SpreadMethod for other algorithms.
"""
function spread(
    points::Vector{CartesianIndex{2}},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real};
    res = 1,
    limit = Inf,
    algorithm = Tomlin(),
)
    spread(algorithm, points, initial, friction; res, limit)
end

function spread(
    points::AbstractMatrix{<:Real},
    initial::AbstractMatrix{<:Real},
    friction::AbstractMatrix{<:Real},
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
