abstract type HorizonMethod end

"""
    GridSweep(directions=8)

Fast sweep method for computing horizon angles.

# Arguments
- `directions::Int`: Number of directions. Must be 4, 8, 16, 32, etc.
  - 4: cardinal (N, E, S, W)
  - 8: cardinal + diagonal
  - 16+: uses rotated grid resampling for intermediate angles
"""
Base.@kwdef struct GridSweep <: HorizonMethod
    directions::Int = 8
end


"""
    horizon_angle(dem; method=GridSweep(), cellsize=cellsize(dem), maxdist=Inf)

Compute horizon angles for each cell in a DEM.

The horizon angle in a given direction is the maximum elevation angle from the cell
to any point along that direction. Positive angles indicate terrain rising above
the horizontal (sheltered), negative angles indicate terrain falling below (exposed).

# Arguments
- `dem`: Digital elevation model matrix
- `method`: `GridSweep(directions=8)` - directions can be 4, 8, 16, 32, ...
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `maxdist`: Maximum search distance in coordinate units

# Returns
3D array of size `(rows, cols, directions)` with horizon angles in degrees.
Direction order: N, NE, E, SE, S, SW, W, NW (for 8 dirs) or N, E, S, W (for 4 dirs).

Angles are in degrees, consistent with `slope` and `aspect`.
"""
function horizon_angle(
    dem::AbstractMatrix{<:Real};
    method::HorizonMethod = GridSweep(),
    cellsize = cellsize(dem),
    maxdist = Inf,
)
    _horizon_angle(method, dem, cellsize, maxdist)
end

function _horizon_angle(method::GridSweep, dem, cellsize, maxdist)
    ndirs = method.directions
    if ndirs == 4
        _horizon_angle_cardinal(dem, cellsize, maxdist)
    elseif ndirs == 8
        _horizon_angle_8(dem, cellsize, maxdist)
    elseif ndirs % 8 == 0 && ndirs > 8
        _horizon_angle_rotated(dem, ndirs, cellsize, maxdist)
    else
        throw(ArgumentError("GridSweep directions must be 4, 8, 16, 32, ..., got $ndirs"))
    end
end

# Single sweep along a line: start at `start`, step by `step`, look back at `-step`
function _sweep_line!(out::AbstractMatrix{Float32}, dem, start::CartesianIndex{2}, step::CartesianIndex{2}, dist::Float64, grid::CartesianIndices{2})
    max_tan = -Inf
    pos = start
    while pos in grid
        elev_pos = Float64(dem[pos])
        if isnan(elev_pos)
            # Reset when crossing NaN - next valid cell is a new edge
            max_tan = -Inf
            pos += step
            continue
        end
        prev = pos - step
        prev_valid = prev in grid && !isnan(Float64(dem[prev]))
        if prev_valid
            tan_angle = (Float64(dem[prev]) - elev_pos) / dist
            max_tan = max(max_tan, tan_angle)
        end
        # If at edge (no valid terrain behind), assume flat horizon (0°)
        out[pos] = Float32(atand(max_tan == -Inf ? 0.0 : max_tan))
        pos += step
    end
end

function _horizon_angle_cardinal(dem, cellsize, maxdist)
    # Returns 3D array with 4 directions: N, E, S, W
    T = Float32
    nrows, ncols = size(dem)
    grid = CartesianIndices(dem)
    result = fill(T(NaN), nrows, ncols, 4)

    cs1, cs2 = Float64(cellsize[1]), Float64(cellsize[2])
    dist_ns = abs(cs1)
    dist_ew = abs(cs2)

    N = view(result, :, :, 1)
    E = view(result, :, :, 2)
    S = view(result, :, :, 3)
    W = view(result, :, :, 4)

    for col in 1:ncols
        _sweep_line!(S, dem, CartesianIndex(nrows, col), CartesianIndex(-1, 0), dist_ns, grid)
        _sweep_line!(N, dem, CartesianIndex(1, col), CartesianIndex(1, 0), dist_ns, grid)
    end
    for row in 1:nrows
        _sweep_line!(E, dem, CartesianIndex(row, ncols), CartesianIndex(0, -1), dist_ew, grid)
        _sweep_line!(W, dem, CartesianIndex(row, 1), CartesianIndex(0, 1), dist_ew, grid)
    end

    return result
end

function _horizon_angle_8(dem, cellsize, maxdist)
    # Returns 3D array with 8 directions: N, NE, E, SE, S, SW, W, NW
    T = Float32
    nrows, ncols = size(dem)
    grid = CartesianIndices(dem)
    result = fill(T(NaN), nrows, ncols, 8)

    cs1, cs2 = Float64(cellsize[1]), Float64(cellsize[2])
    dist_ns = abs(cs1)
    dist_ew = abs(cs2)
    dist_diag = sqrt(cs1^2 + cs2^2)

    N = view(result, :, :, 1)
    NE = view(result, :, :, 2)
    E = view(result, :, :, 3)
    SE = view(result, :, :, 4)
    S = view(result, :, :, 5)
    SW = view(result, :, :, 6)
    W = view(result, :, :, 7)
    NW = view(result, :, :, 8)

    # Sweeps starting from top/bottom edges
    for col in 1:ncols
        _sweep_line!(S, dem, CartesianIndex(nrows, col), CartesianIndex(-1, 0), dist_ns, grid)
        _sweep_line!(N, dem, CartesianIndex(1, col), CartesianIndex(1, 0), dist_ns, grid)
        _sweep_line!(SW, dem, CartesianIndex(nrows, col), CartesianIndex(-1, 1), dist_diag, grid)
        _sweep_line!(NE, dem, CartesianIndex(1, col), CartesianIndex(1, -1), dist_diag, grid)
        _sweep_line!(SE, dem, CartesianIndex(nrows, col), CartesianIndex(-1, -1), dist_diag, grid)
        _sweep_line!(NW, dem, CartesianIndex(1, col), CartesianIndex(1, 1), dist_diag, grid)
    end

    # Sweeps starting from left/right edges
    for row in 1:nrows
        _sweep_line!(E, dem, CartesianIndex(row, ncols), CartesianIndex(0, -1), dist_ew, grid)
        _sweep_line!(W, dem, CartesianIndex(row, 1), CartesianIndex(0, 1), dist_ew, grid)
    end
    for row in 1:nrows-1
        _sweep_line!(SW, dem, CartesianIndex(row, 1), CartesianIndex(-1, 1), dist_diag, grid)
        _sweep_line!(SE, dem, CartesianIndex(row, ncols), CartesianIndex(-1, -1), dist_diag, grid)
        _sweep_line!(NE, dem, CartesianIndex(row + 1, ncols), CartesianIndex(1, -1), dist_diag, grid)
        _sweep_line!(NW, dem, CartesianIndex(row + 1, 1), CartesianIndex(1, 1), dist_diag, grid)
    end

    return result
end

# GPU-compatible bilinear interpolation (no size() calls inside)
function _bilinear_kernel(dem, r, c, nrows, ncols)
    r0 = floor(Int, r)
    c0 = floor(Int, c)
    r1 = min(r0 + 1, nrows)
    c1 = min(c0 + 1, ncols)
    fr = r - r0
    fc = c - c0
    return (1-fr) * (1-fc) * dem[r0, c0] +
           fr * (1-fc) * dem[r1, c0] +
           (1-fr) * fc * dem[r0, c1] +
           fr * fc * dem[r1, c1]
end

function _horizon_angle_rotated(dem, ndirs::Int, cellsize, maxdist)
    nrows, ncols = size(dem)
    cs1, cs2 = Float64(cellsize[1]), Float64(cellsize[2])
    n_rotations = ndirs ÷ 8
    result = zeros(Float32, nrows, ncols, ndirs)

    for rot in 0:n_rotations-1
        angle = rot * (π / 4 / n_rotations)  # Rotation angle in radians

        if rot == 0
            # No rotation needed - use original GridSweep
            h = _horizon_angle_8(dem, cellsize, maxdist)
            for i in 1:8
                di = (i - 1) * n_rotations + 1
                result[:, :, di] .= view(h, :, :, i)
            end
        else
            # Rotate DEM, run GridSweep, map back
            rotated = _rotate_dem(dem, angle, cs1, cs2)
            h = _horizon_angle_8(rotated, cellsize, maxdist)

            # Map results back to original grid positions
            for i in 1:8
                di = (i - 1) * n_rotations + 1 + rot
                _unrotate_result!(view(result, :, :, di), view(h, :, :, i), angle, nrows, ncols)
            end
        end
    end
    return result
end

function _rotate_dem(dem, angle::Float64, cs1, cs2)
    nrows, ncols = size(dem)
    diag = ceil(Int, hypot(nrows, ncols))

    backend = get_backend(dem)
    rotated = KernelAbstractions.allocate(backend, Float64, diag, diag)
    fill!(rotated, NaN)

    cos_a, sin_a = cos(angle), sin(angle)
    center_r, center_c = (nrows + 1) / 2, (ncols + 1) / 2
    center_rot = (diag + 1) / 2

    workgroup = backend isa KernelAbstractions.CPU ? 1 : 256
    kernel! = _rotate_kernel!(backend, workgroup)
    kernel!(rotated, dem, cos_a, sin_a, center_r, center_c, center_rot, nrows, ncols, ndrange=(diag, diag))
    KernelAbstractions.synchronize(backend)
    return rotated
end

@kernel function _rotate_kernel!(rotated, dem, cos_a, sin_a, center_r, center_c, center_rot, nrows, ncols)
    r, c = @index(Global, NTuple)
    dr = r - center_rot
    dc = c - center_rot
    src_r = center_r + dr * cos_a + dc * sin_a
    src_c = center_c - dr * sin_a + dc * cos_a

    if src_r >= 1 && src_r <= nrows && src_c >= 1 && src_c <= ncols
        rotated[r, c] = _bilinear_kernel(dem, src_r, src_c, nrows, ncols)
    end
end

function _unrotate_result!(dest, src, angle::Float64, nrows, ncols)
    diag = size(src, 1)
    cos_a, sin_a = cos(angle), sin(angle)
    center_r, center_c = (nrows + 1) / 2, (ncols + 1) / 2
    center_rot = (diag + 1) / 2

    backend = get_backend(dest)
    workgroup = backend isa KernelAbstractions.CPU ? 1 : 256
    kernel! = _unrotate_kernel!(backend, workgroup)
    kernel!(dest, src, cos_a, sin_a, center_r, center_c, center_rot, diag, ndrange=(nrows, ncols))
    KernelAbstractions.synchronize(backend)
end

@kernel function _unrotate_kernel!(dest, src, cos_a, sin_a, center_r, center_c, center_rot, diag)
    r, c = @index(Global, NTuple)
    dr = r - center_r
    dc = c - center_c
    rot_r = center_rot + dr * cos_a - dc * sin_a
    rot_c = center_rot + dr * sin_a + dc * cos_a

    if rot_r >= 1 && rot_r <= diag && rot_c >= 1 && rot_c <= diag
        val = _bilinear_kernel(src, rot_r, rot_c, diag, diag)
        if isnan(val) || val < -90 || val > 90
            dest[r, c] = 0f0
        else
            dest[r, c] = Float32(val)
        end
    else
        dest[r, c] = 0f0
    end
end

function _bilinear(dem, r::Float64, c::Float64)
    r0 = floor(Int, r)
    c0 = floor(Int, c)
    r1 = min(r0 + 1, size(dem, 1))
    c1 = min(c0 + 1, size(dem, 2))
    fr = r - r0
    fc = c - c0
    return (1-fr) * (1-fc) * dem[r0, c0] +
           fr * (1-fc) * dem[r1, c0] +
           (1-fr) * fc * dem[r0, c1] +
           fr * fc * dem[r1, c1]
end

function _build_max_pyramid(dem)
    nrows, ncols = size(dem)
    nlevels = floor(Int, log2(min(nrows, ncols)))
    levels = Vector{Matrix{Float32}}(undef, nlevels + 1)

    levels[1] = Matrix{Float32}(undef, nrows, ncols)
    for i in eachindex(dem)
        levels[1][i] = Float32(dem[i])
    end

    for lv in 2:nlevels+1
        prev = levels[lv-1]
        nr, nc = size(prev, 1) ÷ 2, size(prev, 2) ÷ 2
        level = Matrix{Float32}(undef, nr, nc)
        for r in 1:nr
            for c in 1:nc
                level[r, c] = max(
                    prev[2r-1, 2c-1], prev[2r-1, 2c],
                    prev[2r, 2c-1], prev[2r, 2c]
                )
            end
        end
        levels[lv] = level
    end
    return levels
end

function _sample_pyramid(pyramid, r::Float64, c::Float64, level::Int)
    lvl = pyramid[level + 1]
    scale = 1 << level
    pr = clamp(ceil(Int, r / scale), 1, size(lvl, 1))
    pc = clamp(ceil(Int, c / scale), 1, size(lvl, 2))
    return lvl[pr, pc]
end

"""
    sky_view_factor(dem; method=GridSweep(), cellsize=cellsize(dem), maxdist=Inf)

Compute the Sky View Factor (SVF) for each cell in a DEM.

SVF is the fraction of the sky hemisphere visible from each point, ranging from 0 (fully
obstructed) to 1 (full sky visible). It is computed as the mean of cos²(horizon_angle)
across all directions.

# Arguments
- `dem`: Digital elevation model matrix
- `method`: Horizon computation method (default: `GridSweep()`)
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `maxdist`: Maximum search distance in coordinate units

# Returns
A matrix of SVF values in the range [0, 1].
"""
function sky_view_factor(
    dem::AbstractMatrix{<:Real};
    method::HorizonMethod = GridSweep(),
    cellsize = cellsize(dem),
    maxdist = Inf,
)
    horizons = horizon_angle(dem; method, cellsize, maxdist)
    ndirs = size(horizons, 3)
    result = zeros(Float32, size(dem))
    for di in 1:ndirs
        for i in CartesianIndices(result)
            result[i] += cosd(horizons[i, di])^2
        end
    end
    result ./= ndirs
    return result
end
