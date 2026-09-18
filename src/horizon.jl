"""
    horizon_angle(dem; directions=16, cellsize=cellsize(dem), observer_height=0.0, missing_elevation=0.0)

Compute horizon angles for each cell in a DEM.

The horizon angle in a given direction is the maximum elevation angle from the cell
to any point along that direction. Positive angles indicate terrain rising above
the horizontal (sheltered), negative angles indicate terrain falling below (exposed).

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `directions`: 16 by default, can be 4, 8, 16, 32, ...
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `observer_height`: Height of the observer's eye above the terrain (default `0.0`). A
  positive value raises the viewpoint, lowering the computed horizon angles.
- `missing_elevation`: Elevation substituted for NaN cells along a sweep ray
  (default `0.0`). Use `Inf` to make NaN fully occluding. Cells whose own
  elevation is NaN always produce `NaN` in the output.

# Returns
3D array of size `(rows, cols, directions)` with horizon angles in degrees.
"""
function horizon_angle(dem::AbstractMatrix{<:Real};
    directions::Int = 16,
    cellsize = cellsize(dem),
    observer_height::Real = 0.0,
    missing_elevation::Real = 0.0,
)
    ispow2(directions) && directions >= 4 ||
        throw(ArgumentError("directions must be 4, 8, 16, 32, ..., got $directions"))
    me = Float64(missing_elevation)
    op = HorizonAngle(Float64(observer_height))
    if directions == 4
        _sweep_directions_cardinal(op, dem, cellsize, me)
    elseif directions == 8
        _sweep_directions_8(op, dem, cellsize, me)
    else
        _sweep_directions_rotated(op, dem, directions, cellsize, me)
    end
end

# Allocate the (rows, cols, directions) output. The default uses `similar` so that
# plain matrices and GPU arrays propagate their storage type. Extensions override
# this for wrapper types like `Raster` and `GeoArray` to attach a direction
# dimension that `similar` would otherwise drop.
_alloc_directions(dem, T, ndirs) =
    similar(dem, T, size(dem, 1), size(dem, 2), ndirs)

# Sweep operations. Singleton structs selected at the call site so the inner
# branch specializes and the unused path is eliminated (no runtime if/else).
abstract type SweepOp end

# Per-cell backward look-back producing the horizon angle (a reduction). `height`
# raises the observer's eye above the terrain (metres).
struct HorizonAngle <: SweepOp
    height::Float64
end

# Per-observer running-max along the line, accumulating the fraction of cells
# each observer can see (used by `total_viewshed`). `height` raises the observer's
# eye above the terrain (metres), so cells on a uniform/grazing slope are seen at a
# progressively shallower depression angle and count as visible.
struct TotalVisibility <: SweepOp
    height::Float64
end
# KernelAbstractions kernels for line sweeps

# Sweep from row edge (top or bottom) - index gives column
@kernel function _sweep_row_edge_kernel!(op::SweepOp, out, @Const(dem), start_row::Int, step_r::Int, step_c::Int, dist::Float64, nrows::Int, ncols::Int, missing_elevation::Float64)
    col = @index(Global)
    _sweep_line_ka!(op, out, dem, start_row, col, step_r, step_c, dist, nrows, ncols, missing_elevation)
end

# Sweep from col edge (left or right) - index gives row
@kernel function _sweep_col_edge_kernel!(op::SweepOp, out, @Const(dem), row_offset::Int, start_col::Int, step_r::Int, step_c::Int, dist::Float64, nrows::Int, ncols::Int, missing_elevation::Float64)
    idx = @index(Global)
    _sweep_line_ka!(op, out, dem, idx + row_offset, start_col, step_r, step_c, dist, nrows, ncols, missing_elevation)
end

# Walk one ray and apply `op` per cell, sharing the same backward look-back.
#
# For each cell, scan every prior cell on the ray (nearest first) tracking the
# running max of `Δz / (k·dist)`. NaN look-back substitutes `missing_elevation`;
# NaN at the current cell leaves the output as NaN (from the initial fill). A
# prior cell is *visible* exactly when its angle sets a new running max, so:
#
# - `HorizonAngle` writes `atand(max_tan)` (the horizon angle).
# - `TotalVisibility` counts those visible cells and writes the fraction of the
#   cells behind it that are visible.
#
# The `op isa ...` tests are resolved at compile time, so each specialization
# keeps a minimal hot loop.
function _sweep_line_ka!(op::SweepOp, out, dem, start_r, start_c, step_r, step_c, dist, nrows, ncols, missing_elevation::Float64)
    r, c = start_r, start_c
    step_count = 0
    @inbounds while r >= 1 && r <= nrows && c >= 1 && c <= ncols
        elev_pos = Float64(dem[r, c])
        if !isnan(elev_pos)
            # Raise the observer's eye above the terrain (0 reproduces the
            # ground-level horizon/visibility).
            eye = elev_pos + op.height
            max_tan = -Inf
            visible = 0
            valid = 0
            for k in 1:step_count
                pr = r - k * step_r
                pc = c - k * step_c
                prev_elev = Float64(dem[pr, pc])
                is_missing = isnan(prev_elev)
                if is_missing
                    prev_elev = missing_elevation
                end
                tan_angle = (prev_elev - eye) / (k * dist)
                # Count only real cells as visibility targets, but still let
                # missing cells feed the running max (so `Inf` occludes), exactly
                # as a missing cell is ignored by the horizon max.
                if op isa TotalVisibility && !is_missing
                    valid += 1
                    if tan_angle >= max_tan
                        visible += 1
                    end
                end
                if tan_angle > max_tan
                    max_tan = tan_angle
                end
            end
            if op isa TotalVisibility
                out[r, c] = valid == 0 ? 1.0f0 : Float32(visible / valid)
            else
                out[r, c] = Float32(atand(max_tan == -Inf ? 0.0 : max_tan))
            end
        end
        r += step_r
        c += step_c
        step_count += 1
    end
end

function _sweep_directions_cardinal(op::SweepOp, dem, cellsize, missing_elevation::Float64)
    # Returns 3D array with 4 directions: N, E, S, W
    T = Float32
    nrows, ncols = size(dem)
    result = _alloc_directions(dem, T, 4)
    fill!(result, T(NaN))
    backend = get_backend(parent(dem))

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
    sweep_row_edge!(op, S, dem, nrows, -1, 0, dist_ns, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, N, dem, 1, 1, 0, dist_ns, nrows, ncols, missing_elevation, ndrange=ncols)
    # E/W from left/right edges (step_r=0 for cardinal)
    sweep_col_edge!(op, E, dem, 0, ncols, 0, -1, dist_ew, nrows, ncols, missing_elevation, ndrange=nrows)
    sweep_col_edge!(op, W, dem, 0, 1, 0, 1, dist_ew, nrows, ncols, missing_elevation, ndrange=nrows)

    KernelAbstractions.synchronize(backend)
    return result
end

function _sweep_directions_8(op::SweepOp, dem, cellsize, missing_elevation::Float64)
    # Returns 3D array with 8 directions: N, NE, E, SE, S, SW, W, NW
    T = Float32
    nrows, ncols = size(dem)
    result = _alloc_directions(dem, T, 8)
    fill!(result, T(NaN))
    backend = get_backend(parent(dem))

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
    sweep_row_edge!(op, S, dem, nrows, -1, 0, dist_ns, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, N, dem, 1, 1, 0, dist_ns, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, SW, dem, nrows, -1, 1, dist_diag, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, NE, dem, 1, 1, -1, dist_diag, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, SE, dem, nrows, -1, -1, dist_diag, nrows, ncols, missing_elevation, ndrange=ncols)
    sweep_row_edge!(op, NW, dem, 1, 1, 1, dist_diag, nrows, ncols, missing_elevation, ndrange=ncols)

    # Sweeps from left/right edges (ndrange=nrows for E/W)
    sweep_col_edge!(op, E, dem, 0, ncols, 0, -1, dist_ew, nrows, ncols, missing_elevation, ndrange=nrows)
    sweep_col_edge!(op, W, dem, 0, 1, 0, 1, dist_ew, nrows, ncols, missing_elevation, ndrange=nrows)

    # Diagonal sweeps from left/right edges (ndrange=nrows-1)
    sweep_col_edge!(op, SW, dem, 0, 1, -1, 1, dist_diag, nrows, ncols, missing_elevation, ndrange=nrows-1)
    sweep_col_edge!(op, SE, dem, 0, ncols, -1, -1, dist_diag, nrows, ncols, missing_elevation, ndrange=nrows-1)
    sweep_col_edge!(op, NE, dem, 1, ncols, 1, -1, dist_diag, nrows, ncols, missing_elevation, ndrange=nrows-1)
    sweep_col_edge!(op, NW, dem, 1, 1, 1, 1, dist_diag, nrows, ncols, missing_elevation, ndrange=nrows-1)

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

function _sweep_directions_rotated(op::SweepOp, dem, ndirs::Int, cellsize, missing_elevation::Float64)
    nrows, ncols = size(dem)
    cs1, cs2 = Float64(cellsize[1]), Float64(cellsize[2])
    n_rotations = ndirs ÷ 8
    result = _alloc_directions(dem, Float32, ndirs)
    fill!(result, 0f0)

    # No rotation needed for first 8
    h = _sweep_directions_8(op, dem, cellsize, missing_elevation)
    for i in 1:8
        di = (i - 1) * n_rotations + 1
        copyto!(view(result, :, :, di), view(h, :, :, i))
    end

    # For the rest we rotate DEM, run, map back
    for rot in 1:n_rotations-1
        angle = rot * (π / 4 / n_rotations)  # Rotation angle in radians

        rotated = _rotate_dem(dem, angle, cs1, cs2)
        h = _sweep_directions_8(op, rotated, cellsize, missing_elevation)

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

    rotated = similar(parent(dem), Float64, diag, diag)
    backend = get_backend(rotated)
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
    sky_view_factor(dem; directions=16, cellsize=cellsize(dem), observer_height=0.0, missing_elevation=0.0)

Compute the Sky View Factor (SVF) for each cell in a DEM.

SVF is the fraction of the sky hemisphere visible from each point, ranging from 0 (fully
obstructed) to 1 (full sky visible). It is computed as the mean of cos²(horizon_angle)
across all directions.

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `directions`: Number of directions (default: 16)
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `observer_height`: Height of the observer's eye above the terrain (default `0.0`),
  forwarded to [`horizon_angle`](@ref).
- `missing_elevation`: Elevation substituted for missing cells in `dem` when computing
  horizons (default `0.0`). See [`horizon_angle`](@ref).

# Returns
A matrix of SVF values in the range [0, 1].
"""
function sky_view_factor(dem::AbstractMatrix{<:Real};
    directions::Int = 16,
    cellsize = cellsize(dem),
    observer_height::Real = 0.0,
    missing_elevation::Real = 0.0,
)
    horizons = horizon_angle(dem; directions, cellsize, observer_height, missing_elevation)
    ndirs = size(horizons, 3)
    result = similar(dem, Float32)
    backend = get_backend(parent(horizons))

    workgroup = backend isa KernelAbstractions.CPU ? 1 : (16, 16)
    kernel! = _sky_view_factor_kernel!(backend, workgroup)
    kernel!(result, horizons, ndirs, ndrange=size(dem))
    KernelAbstractions.synchronize(backend)

    return result
end

@kernel function _total_viewshed_kernel!(result, @Const(vis), ndirs::Int)
    i, j = @index(Global, NTuple)
    res = 0.0f0
    @inbounds for di in 1:ndirs
        res += vis[i, j, di]
    end
    result[i, j] = res / ndirs
end

"""
    total_viewshed(dem; directions=16, cellsize=cellsize(dem), observer_height=0.0, missing_elevation=0.0)

Compute the total viewshed (a normalized visibility index) for each cell in a DEM.

For every cell, and along each direction, the fraction of the cells behind it (toward the
sweep origin) that are directly visible (unobstructed line of sight) is computed; these
fractions are then averaged over all directions. Values range from 0 (nothing visible) to 1
(everything visible). Flat terrain yields 1.0 everywhere.

# Arguments
- `dem`: Digital elevation model matrix

# Keywords
- `directions`: Number of directions (default: 16), can be 4, 8, 16, 32, ...
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `observer_height`: Height of the observer's eye above the terrain (default `0.0`). A small
  positive value (e.g. `1.6` m) lets the observer see over flat/uniform terrain and makes
  the result robust to grazing-angle round-off.
- `missing_elevation`: Elevation substituted for NaN cells along a sweep ray (default `0.0`).
  Missing cells still occlude (feed the running maximum) but are not counted as visibility
  targets. See [`horizon_angle`](@ref).

# Returns
A matrix of total-viewshed values in the range [0, 1].
"""
function total_viewshed(dem::AbstractMatrix{<:Real};
    directions::Int = 16,
    cellsize = cellsize(dem),
    observer_height::Real = 0.0,
    missing_elevation::Real = 0.0,
)
    ispow2(directions) && directions >= 4 ||
        throw(ArgumentError("directions must be 4, 8, 16, 32, ..., got $directions"))
    me = Float64(missing_elevation)
    op = TotalVisibility(Float64(observer_height))
    vis = if directions == 4
        _sweep_directions_cardinal(op, dem, cellsize, me)
    elseif directions == 8
        _sweep_directions_8(op, dem, cellsize, me)
    else
        _sweep_directions_rotated(op, dem, directions, cellsize, me)
    end

    ndirs = size(vis, 3)
    result = similar(dem, Float32)
    backend = get_backend(parent(vis))

    workgroup = backend isa KernelAbstractions.CPU ? 1 : (16, 16)
    kernel! = _total_viewshed_kernel!(backend, workgroup)
    kernel!(result, vis, ndirs, ndrange=size(dem))
    KernelAbstractions.synchronize(backend)

    return result
end

"""
    viewshed(dem, observer; cellsize=cellsize(dem), observer_height=0.0, missing_elevation=0.0)

Compute the viewshed of a single `observer` cell in a DEM.

For every cell, the line of sight from the observer to that cell is traced. Along the line
the terrain elevation is sampled at the points where the sightline crosses each grid
row/column, using bilinear interpolation. The target cell is visible if its elevation angle
from the observer is at least the running maximum of those intermediate angles. This is the
same running-maximum test that `horizon_angle` uses, read as a boolean state. The observer
cell itself is always visible.

# Arguments
- `dem`: Digital elevation model matrix
- `observer`: Observer cell, e.g. a `CartesianIndex` or `(row, col)` tuple

# Keywords
- `cellsize`: Cell size as `(row_size, col_size)` tuple
- `observer_height`: Height of the observer's eye above the terrain (default `0.0`). A small
  positive value (e.g. `1.6` m) lets the observer see over flat ground and grazing slopes.
- `missing_elevation`: Elevation substituted for NaN cells along a line of sight
  (default `0.0`). See [`horizon_angle`](@ref).

# Returns
A `Bool` matrix that is `true` where the cell is visible from the observer. Cells whose own
elevation is NaN (and all cells when the observer elevation is NaN) are `false`.
"""
function viewshed(dem::AbstractMatrix{<:Real}, observer;
    cellsize = cellsize(dem),
    observer_height::Real = 0.0,
    missing_elevation::Real = 0.0,
)
    nrows, ncols = size(dem)
    or, oc = Int(observer[1]), Int(observer[2])
    (1 <= or <= nrows && 1 <= oc <= ncols) ||
        throw(ArgumentError("observer ($or, $oc) is outside the DEM"))
    cs1, cs2 = abs(Float64(cellsize[1])), abs(Float64(cellsize[2]))
    me = Float64(missing_elevation)
    oh = Float64(observer_height)

    result = similar(dem, Bool)
    backend = get_backend(parent(dem))

    workgroup = backend isa KernelAbstractions.CPU ? 1 : (16, 16)
    kernel! = _viewshed_kernel!(backend, workgroup)
    kernel!(result, dem, or, oc, cs1, cs2, nrows, ncols, me, oh, ndrange=size(dem))
    KernelAbstractions.synchronize(backend)

    return result
end

@kernel function _viewshed_kernel!(out, @Const(dem), or::Int, oc::Int, cs1::Float64, cs2::Float64, nrows::Int, ncols::Int, missing_elevation::Float64, observer_height::Float64)
    i, j = @index(Global, NTuple)
    out[i, j] = _line_of_sight(dem, or, oc, i, j, cs1, cs2, nrows, ncols, missing_elevation, observer_height)
end

# Trace the line of sight observer -> (i, j). Sample terrain at the sightline's
# crossing point on each intervening grid line by bilinear interpolation, take the
# running max of the intermediate elevation angles, and report visibility of the
# target. The observer's eye is raised by `observer_height`. NaN samples substitute
# `missing_elevation`; a NaN observer or target is not visible.
@inline function _line_of_sight(dem, or, oc, i, j, cs1, cs2, nrows, ncols, missing_elevation::Float64, observer_height::Float64)
    @inbounds z_obs = Float64(dem[or, oc])
    @inbounds z_target = Float64(dem[i, j])
    (isnan(z_obs) || isnan(z_target)) && return false
    z_eye = z_obs + observer_height

    dr = i - or
    dc = j - oc
    (dr == 0 && dc == 0) && return true

    nsteps = max(abs(dr), abs(dc))
    D = sqrt((dr * cs1)^2 + (dc * cs2)^2)

    run_max = -Inf
    @inbounds for k in 1:nsteps-1
        t = k / nsteps
        z = _bilinear_kernel(dem, or + t * dr, oc + t * dc, nrows, ncols)
        if isnan(z)
            z = missing_elevation
        end
        ang = (z - z_eye) / (t * D)
        if ang > run_max
            run_max = ang
        end
    end

    target_ang = (z_target - z_eye) / D
    return target_ang >= run_max
end
