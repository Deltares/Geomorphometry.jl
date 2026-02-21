"""
    horizon_angle(dem; method=GridSweep(), cellsize=cellsize(dem), maxdist=Inf)

Compute horizon angles for each cell in a DEM.

The horizon angle in a given direction is the maximum elevation angle from the cell
to any point along that direction. Positive angles indicate terrain rising above
the horizontal (sheltered), negative angles indicate terrain falling below (exposed).

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `directions`: 16 by default, can be 4, 8, 16, 32, ...
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `maxdist`: Maximum search distance in coordinate units

# Returns
3D array of size `(rows, cols, directions)` with horizon angles in degrees.
Direction order: N, NE, E, SE, S, SW, W, NW (for 8 dirs) or N, E, S, W (for 4 dirs).

Angles are in degrees, consistent with `slope` and `aspect`.
"""
function horizon_angle(
    dem::AbstractMatrix{<:Real};
    directions::Int = 16,
    cellsize = cellsize(dem),
    maxdist = Inf,
)
    if directions == 4
        _horizon_angle_cardinal(dem, cellsize, maxdist)
    elseif directions == 8
        _horizon_angle_8(dem, cellsize, maxdist)
    elseif directions % 8 == 0 && directions > 8
        _horizon_angle_rotated(dem, directions, cellsize, maxdist)
    else
        throw(ArgumentError("directions must be 4, 8, 16, 32, ..., got $directions"))
    end
end

# KernelAbstractions kernels for line sweeps

# Sweep from row edge (top or bottom) - index gives column
@kernel function _sweep_row_edge_kernel!(out, @Const(dem), start_row::Int, step_r::Int, step_c::Int, dist::Float64, nrows::Int, ncols::Int)
    col = @index(Global)
    _sweep_line_ka!(out, dem, start_row, col, step_r, step_c, dist, nrows, ncols)
end

# Sweep from col edge (left or right) - index gives row
@kernel function _sweep_col_edge_kernel!(out, @Const(dem), row_offset::Int, start_col::Int, step_r::Int, step_c::Int, dist::Float64, nrows::Int, ncols::Int)
    idx = @index(Global)
    _sweep_line_ka!(out, dem, idx + row_offset, start_col, step_r, step_c, dist, nrows, ncols)
end

# Shared sweep implementation for any direction
function _sweep_line_ka!(out, dem, start_r, start_c, step_r, step_c, dist, nrows, ncols)
    max_tan = -Inf
    r, c = start_r, start_c
    @inbounds while r >= 1 && r <= nrows && c >= 1 && c <= ncols
        elev_pos = Float64(dem[r, c])
        if isnan(elev_pos)
            max_tan = -Inf
        else
            pr, pc = r - step_r, c - step_c
            if pr >= 1 && pr <= nrows && pc >= 1 && pc <= ncols
                prev_elev = Float64(dem[pr, pc])
                if !isnan(prev_elev)
                    tan_angle = (prev_elev - elev_pos) / dist
                    max_tan = max(max_tan, tan_angle)
                end
            end
            out[r, c] = Float32(atand(max_tan == -Inf ? 0.0 : max_tan))
        end
        r += step_r
        c += step_c
    end
end

function _horizon_angle_cardinal(dem, cellsize, maxdist)
    # Returns 3D array with 4 directions: N, E, S, W
    T = Float32
    nrows, ncols = size(dem)
    backend = get_backend(dem)
    result = KernelAbstractions.allocate(backend, T, nrows, ncols, 4)
    fill!(result, T(NaN))

    cs1, cs2 = Float64(cellsize[1]), Float64(cellsize[2])
    dist_ns = abs(cs1)
    dist_ew = abs(cs2)

    N = view(result, :, :, 1)
    E = view(result, :, :, 2)
    S = view(result, :, :, 3)
    W = view(result, :, :, 4)

    workgroup = backend isa KernelAbstractions.CPU ? 1 : 256
    sweep_row_edge! = _sweep_row_edge_kernel!(backend, workgroup)
    sweep_col_edge! = _sweep_col_edge_kernel!(backend, workgroup)

    # N/S from top/bottom edges (step_c=0 for cardinal)
    sweep_row_edge!(S, dem, nrows, -1, 0, dist_ns, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(N, dem, 1, 1, 0, dist_ns, nrows, ncols, ndrange=ncols)
    # E/W from left/right edges (step_r=0 for cardinal)
    sweep_col_edge!(E, dem, 0, ncols, 0, -1, dist_ew, nrows, ncols, ndrange=nrows)
    sweep_col_edge!(W, dem, 0, 1, 0, 1, dist_ew, nrows, ncols, ndrange=nrows)

    KernelAbstractions.synchronize(backend)
    return result
end

function _horizon_angle_8(dem, cellsize, maxdist)
    # Returns 3D array with 8 directions: N, NE, E, SE, S, SW, W, NW
    T = Float32
    nrows, ncols = size(dem)
    backend = get_backend(dem)
    result = KernelAbstractions.allocate(backend, T, nrows, ncols, 8)
    fill!(result, T(NaN))

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

    workgroup = backend isa KernelAbstractions.CPU ? 1 : 256
    sweep_row_edge! = _sweep_row_edge_kernel!(backend, workgroup)
    sweep_col_edge! = _sweep_col_edge_kernel!(backend, workgroup)

    # Sweeps from top/bottom edges (ndrange=ncols)
    sweep_row_edge!(S, dem, nrows, -1, 0, dist_ns, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(N, dem, 1, 1, 0, dist_ns, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(SW, dem, nrows, -1, 1, dist_diag, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(NE, dem, 1, 1, -1, dist_diag, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(SE, dem, nrows, -1, -1, dist_diag, nrows, ncols, ndrange=ncols)
    sweep_row_edge!(NW, dem, 1, 1, 1, dist_diag, nrows, ncols, ndrange=ncols)

    # Sweeps from left/right edges (ndrange=nrows for E/W)
    sweep_col_edge!(E, dem, 0, ncols, 0, -1, dist_ew, nrows, ncols, ndrange=nrows)
    sweep_col_edge!(W, dem, 0, 1, 0, 1, dist_ew, nrows, ncols, ndrange=nrows)

    # Diagonal sweeps from left/right edges (ndrange=nrows-1)
    sweep_col_edge!(SW, dem, 0, 1, -1, 1, dist_diag, nrows, ncols, ndrange=nrows-1)
    sweep_col_edge!(SE, dem, 0, ncols, -1, -1, dist_diag, nrows, ncols, ndrange=nrows-1)
    sweep_col_edge!(NE, dem, 1, ncols, 1, -1, dist_diag, nrows, ncols, ndrange=nrows-1)
    sweep_col_edge!(NW, dem, 1, 1, 1, 1, dist_diag, nrows, ncols, ndrange=nrows-1)

    KernelAbstractions.synchronize(backend)
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
    backend = get_backend(dem)
    result = KernelAbstractions.allocate(backend, Float32, nrows, ncols, ndirs)
    fill!(result, 0f0)

    # No rotation needed for first 8
    h = _horizon_angle_8(dem, cellsize, maxdist)
    for i in 1:8
        di = (i - 1) * n_rotations + 1
        copyto!(view(result, :, :, di), view(h, :, :, i))
    end

    # For the rest we rotate DEM, run, map back
    for rot in 1:n_rotations-1
        angle = rot * (π / 4 / n_rotations)  # Rotation angle in radians

        rotated = _rotate_dem(dem, angle, cs1, cs2)
        h = _horizon_angle_8(rotated, cellsize, maxdist)

        # Map results back to original grid positions
        for i in 1:8
            di = (i - 1) * n_rotations + 1 + rot
            _unrotate_result!(view(result, :, :, di), view(h, :, :, i), angle, nrows, ncols)
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

@kernel function _sky_view_factor_kernel!(result, @Const(horizons), ndirs::Int)
    i, j = @index(Global, NTuple)
    res = 0.0f0
    @inbounds for di in 1:ndirs
        res += cosd(horizons[i, j, di])^2
    end
    result[i, j] = res / ndirs
end

"""
    sky_view_factor(dem; directions=16, cellsize=cellsize(dem), maxdist=Inf)

Compute the Sky View Factor (SVF) for each cell in a DEM.

SVF is the fraction of the sky hemisphere visible from each point, ranging from 0 (fully
obstructed) to 1 (full sky visible). It is computed as the mean of cos²(horizon_angle)
across all directions.

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `directions`: Number of directions (default: 16)
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `maxdist`: Maximum search distance in coordinate units

# Returns
A matrix of SVF values in the range [0, 1].
"""
function sky_view_factor(
    dem::AbstractMatrix{<:Real};
    directions::Int = 16,
    cellsize = cellsize(dem),
    maxdist = Inf,
)
    horizons = horizon_angle(dem; directions, cellsize, maxdist)
    ndirs = size(horizons, 3)
    backend = get_backend(horizons)
    result = KernelAbstractions.allocate(backend, Float32, size(dem))

    workgroup = backend isa KernelAbstractions.CPU ? 1 : (16, 16)
    kernel! = _sky_view_factor_kernel!(backend, workgroup)
    kernel!(result, horizons, ndirs, ndrange=size(dem))
    KernelAbstractions.synchronize(backend)

    return result
end
