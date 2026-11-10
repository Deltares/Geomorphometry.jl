"""
    eachcell(dem)

Iterate over every cell of `dem`, as indices that `dem` accepts directly.
Matrices yield `CartesianIndex{2}`, other grids fall back to `eachindex`.
Grids with their own cell type extend this method, together with
[`neighbors`](@ref).
"""
eachcell(dem) = eachindex(dem)
eachcell(dem::AbstractMatrix) = CartesianIndices(dem)

"""
    neighbors(dem)

Iterate over `(cell, neighbors(dem, cell))` for every cell in `dem`. Neighbor
iterators contain only indices within `dem`.
"""
neighbors(dem) = ((cell, neighbors(dem, cell)) for cell in eachcell(dem))

_neighborvalues(dem, cell) = (dem[neighbor] for neighbor in neighbors(dem, cell))

"""
    mapneighbors!(f, dst, dem; order = nothing, threaded = false)

Apply `f(cell, value, neighbor_values)` to every cell in `dem`, writing results
directly into `dst` and returning it. `neighbor_values` contains samples from
[`neighbors`](@ref). Calls to `f` must be independent: backends may traverse in
`order` or run calls concurrently when `threaded` is `true`. The fallback
ignores both keywords.
"""
function mapneighbors!(f::F, dst, dem; order=nothing, threaded=false) where {F}
    for cell in eachcell(dem)
        dst[cell] = f(cell, dem[cell], _neighborvalues(dem, cell))
    end
    return dst
end

"""
    mapneighbors(f, dem; order = nothing, threaded = false)

Apply `f(cell, value, neighbor_values)` to every cell in `dem` and return the
results in a grid shaped like `dem`. `neighbor_values` contains samples from
[`neighbors`](@ref). Calls to `f` must be independent: backends may traverse in
`order` or run calls concurrently when `threaded` is `true`. The fallback
ignores both keywords.
"""
function mapneighbors(f::F, dem; order=nothing, threaded=false) where {F}
    cells = eachcell(dem)
    values = map(cell -> f(cell, dem[cell], _neighborvalues(dem, cell)), cells)
    dst = similar(dem, eltype(values))
    for (cell, value) in zip(cells, values)
        dst[cell] = value
    end
    return dst
end
