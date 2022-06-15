const sqrt2 = sqrt(2.0)
const distance_8 = @SMatrix[sqrt2 1 sqrt2; 1 Inf 1; sqrt2 1 sqrt2]
const distance_4 = @SMatrix[Inf 1 Inf; 1 Inf 1; Inf 1 Inf]
const Δ = CartesianIndex(1, 1)

"""
```
spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf)
```
Total friction distance spread from `points` as by *Tomlin (1983)* [^tomlin1983].
This is also the method implemented by [PCRaster](https://pcraster.geo.uu.nl/pcraster/4.0.2/doc/manual/op_spread.html).

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Matrix{<:Real}` Input Array
- `initial::Matrix{<:Real}` Factor to exaggerate elevation
- `friction::Matrix{<:Real}` Resolution of cell size
- `res=1` Resolution or cell size
- `limit=Inf` Initial fill value

[^tomlin1983]: Tomlin, Charles Dana. 1983. Digital Cartographic Modeling Techniques in Environmental Planning. Yale University.
"""
function spread(points::AbstractMatrix{<:Real}, initial::AbstractMatrix{<:Real}, friction::AbstractMatrix{<:Real}; res=1, limit=Inf)

    ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    ofriction[begin+1:end-1, begin+1:end-1] .= friction

    result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    r = @view result[1:end-1, 1:end-1]
    locations = points .> 0
    r[locations] .= initial[locations]

    # Construct stack for locations
    mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    mask[begin+1:end-1, begin+1:end-1] .= false

    II = CartesianIndices(size(points))
    stack = Deque{CartesianIndex}()
    for I in II[locations]
        push!(stack, I)
    end

    # Step 1: Set the distance of the starting node to 0 and the distances of all other nodes to the highest value possible.
    sdata = zeros(MMatrix{3,3})
    mcell = MMatrix{3,3}(false, false, false, false, false, false, false, false, false)

    # Step 3: For each of the active node’s adjacent neighbors, set its distance to whichever is
    # less: its current distance value or the sum of the distance of the active node plus the
    # weight of the arc from the active node to that neighbor.
    while !isempty(stack)
        spread!(stack, mask, result, ofriction, sdata, mcell, res)
    end
    r
end


function spread!(stack, mask, result, ofriction, sdata, mcell, res)
    I = popfirst!(stack)
    mask[I] = true
    patch = I-Δ:I+Δ

    rdata = view(result, patch)
    fdata = view(ofriction, patch)

    # New distance is cell_distance + average friction values
    for i ∈ eachindex(sdata)
        sdata[i] = muladd(fdata[i] + fdata[2, 2], res / 2 * distance_8[i], rdata[2, 2])
        mcell[i] = sdata[i] < rdata[i]  # cells where new distance is lower
    end
    rdata[mcell] .= sdata[mcell]
    result[patch] .= rdata

    # Add new cells to stack
    for I in patch[mcell]
        mask[I] || push!(stack, I)
    end
end

"""
```
spread2(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Matrix{<:Real}; res=1, limit=Inf, iterations=3)
```

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
- `iterations=3` Number of pushbroom iterations

[^eastman1989]: Eastman, J. Ronald. 1989. ‘Pushbroom Algorithms for Calculating Distances in Raster Grids’. In Proceedings, Autocarto, 9:288–97.
"""
function spread2(points::AbstractMatrix{<:Real}, initial::AbstractMatrix{<:Real}, friction::AbstractMatrix{<:Real}; res=1, limit=Inf, iterations=3)

    ofriction = OffsetMatrix(fill(Inf, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    ofriction[begin+1:end-1, begin+1:end-1] .= friction

    result = OffsetMatrix(fill(limit, size(friction) .+ 2), UnitRange.(0, size(points) .+ 1))
    r = @view result[1:end-1, 1:end-1]
    locations = points .> 0
    r[locations] .= initial[locations]

    mask = OffsetMatrix(trues(size(points) .+ 2), UnitRange.(0, size(points) .+ 1))
    mask[begin+1:end-1, begin+1:end-1] .= false

    minval, minidx = [0.0], [CartesianIndex(1, 1)]
    x = @MMatrix zeros(3, 3)

    indices = CartesianIndices(size(points))
    for i in 1:iterations
        II = (i % 2 == 1) ? indices : reverse(indices)
        for I ∈ II
            patch = I-Δ:I+Δ

            rdata = view(result, patch)
            fdata = view(ofriction, patch)

            x .= (fdata .+ fdata[2, 2]) .* res ./ 2 .* distance_8 .+ rdata
            findmin!(minval, minidx, x)
            rdata[2, 2] = min(rdata[2, 2], minval[1])
        end
    end
    r
end

"""
```
spread(points::Matrix{<:Real}, initial::Real, friction::Real; distance=Euclidean(), res=1.0)
```
"""
function spread(points::AbstractMatrix{<:Real}, initial::Real, friction::Real; distance=Euclidean(), res=1.0)
    init = fill(initial, size(points))
    return spread(points, init, friction; distance, res)
end


"""
```
spread(points::Matrix{<:Real}, initial::Matrix{<:Real}, friction::Real; distance=Euclidean(), res=1.0)
```
Optimized (and more accurate) function based on the same friction everywhere.

When the friction is the same everywhere, there's no need for searching the shortest cost path,
as one can just take a direct line to the input points.

The calculated cost is more accurate, as there's no 'zigzag' from cell center to cell center.
"""
function spread(points::AbstractMatrix{<:Real}, initial::AbstractMatrix{<:Real}, friction::Real; distance=Euclidean(), res=1.0)
    locations = points .> 0
    I = CartesianIndices(size(points))

    result = fill(Inf, size(points))
    for location ∈ I[locations]
        for cell ∈ I
            result[cell] = min(evaluate(distance, location.I, cell.I) * res * friction + initial[location], result[cell])
        end
    end
    m = .~isfinite.(points)
    result[m] .= points[m]
    return result
end
