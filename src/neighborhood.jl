# Fallbacks visit every index and delegate lookup to `neighbors(dem, cell)`.

_neighborhood_cells(dem) = eachindex(dem)
_neighborhood_cells(dem::AbstractMatrix) = CartesianIndices(dem)

"""
    neighbors(dem)

Iterate over `(cell, neighbors(dem, cell))` for every cell in `dem`. Neighbor
iterators contain only indices within `dem`.
"""
neighbors(dem) = ((cell, neighbors(dem, cell)) for cell in _neighborhood_cells(dem))

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
    for cell in _neighborhood_cells(dem)
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
    cells = _neighborhood_cells(dem)
    values = map(cell -> f(cell, dem[cell], _neighborvalues(dem, cell)), cells)
    dst = similar(dem, eltype(values))
    for (cell, value) in zip(cells, values)
        dst[cell] = value
    end
    return dst
end
