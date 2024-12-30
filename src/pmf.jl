"""
    B, flags = pmf(A; ωₘ, slope, dhₘ, dh₀, cellsize, adjust, erode)

Applies the progressive morphological filter by *Zhang et al. (2003)* [^keqizhangProgressiveMorphologicalFilter2003](@cite) to `A`.

# Output
- `B::Array{T,2}` Maximum allowable values
- `flags::Array{Float64,2}` A sized array with window sizes if filtered, zero if not filtered.

Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.

# Arguments
- `A::Array{T,2}` Input Array
- `ωₘ::Real=20.` Maximum window size [m]
- `slope::Real=0.01` Terrain slope [m/m]
- `dhₘ::Real=2.5` Maximum elevation threshold [m]
- `dh₀::Real=0.2` Initial elevation threshold [m]
- `cellsize::Real=1.` Cellsize in [m]

# [^zhang2003]: Zhang, Keqi, Shu-Ching Chen, Dean Whitman, Mei-Ling Shyu, Jianhua Yan, and Chengcui Zhang. “A Progressive Morphological Filter for Removing Nonground Measurements from Airborne LIDAR Data.” IEEE Transactions on Geoscience and Remote Sensing 41, no. 4 (2003): 872–82. <https://doi.org/10.1109/TGRS.2003.810682>.
"""
function pmf(
    A::AbstractMatrix{<:Real};
    ωₘ = 20.0,
    slope = 0.01,
    dhₘ = 2.5,
    dh₀ = 0.2,
    cellsize = 1.0,
    circular = false,
    adjust = false,
    erode = false,
)
    _pmf(A, ωₘ, slope, dhₘ, dh₀, cellsize, circular, adjust, erode)
end

function _pmf(
    A::AbstractMatrix{<:Real},
    ωₘ::Real,
    slope::Real,
    dhₘ::Real,
    dh₀::Real,
    cellsize::Real,
    circular::Bool,
    adjust::Bool,
    erode::Bool,
)
    _pmf(
        A,
        ωₘ,
        Fill(slope, size(A)),
        dhₘ,
        Fill(dh₀, size(A)),
        cellsize,
        circular,
        adjust,
        erode,
    )
end

function _pmf(
    A::AbstractMatrix{<:Real},
    ωₘ::Real,
    slope::AbstractMatrix{<:Real},
    dhₘ::Real,
    dh₀::AbstractMatrix{<:Real},
    cellsize::Real,
    circular::Bool,
    adjust::Bool,
    erode::Bool,
)

    # Compute windowsizes and thresholds
    ωₘ = round(Int, ωₘ / cellsize)
    κ_max = floor(Int, log2(ωₘ - 1))  # determine # iterations based on exp growth
    windowsizes = Int.(exp2.(1:κ_max)) .+ 1
    @info windowsizes

    # Compute tresholds
    dwindows = vcat(windowsizes[1], windowsizes)  # prepend first element so we get 0 as diff
    window_diffs = [dwindows[i] - dwindows[i - 1] for i in 2:length(dwindows)]
    # height_tresholds = [min(dhₘ, slope * window_diff * cellsize + dh₀) for window_diff in window_diffs]

    # Set up arrays
    Af = copy(A)  # array to be morphed
    nan_mask = isnan.(Af)
    Af[nan_mask] .= Inf  # Replace NaN with Inf, as to always filter these

    B = copy(A)  # max_elevation raster
    out = copy(A)  # max_elevation raster

    flags = similar(A, Float64)  # 0 = ground, other values indicate window size
    fill!(flags, 0.0)
    flags[nan_mask] .= NaN

    mask = falses(size(A))

    # Iterate over window sizes and height tresholds
    for (i, ωₖ) in enumerate(windowsizes)
        s = (i > 1) && adjust ? dilate(slope, window_diffs[i]) : slope
        @debug "Window $(ωₖ), $(window_diffs[i]) slope sum: $(sum(s))"
        dhₜ = min.(dhₘ, s * window_diffs[i] * cellsize .+ dh₀)
        if erode
            if circular
                mapwindowcirc_approx!(minimum_mask, A, ωₖ, Af, Inf)
            else
                # mapwindow_stack!(minimum, A, ωₖ, Af)
                LocalFilters.erode!(Af, A, ωₖ)
            end
        else
            if circular
                opening_circ!(Af, ωₖ, out)
            else
                LocalFilters.opening!(Af, out, Af, ωₖ)
            end
        end
        mask .= (A .- Af) .> dhₜ
        @info sum(mask)
        for I in eachindex(flags)
            if mask[I] && (flags[I] <= 0)
                flags[I] = ωₖ
            end
        end
        B .= min.(B, Af .+ dhₜ)
    end

    B, flags
end

function pmf(A::AbstractArray{<:Real, 3}; kwargs...)
    size(A, 3) == 1 || throw(ArgumentError("Only singleton 3rd dimension allowed"))
    pmf(view(A, :, :, 1); kwargs...)
end

function pmf2(
    A::AbstractMatrix{<:Real};
    ωₘ = 20.0,
    slope = 0.01,
    dhₘ = 2.5,
    dh₀ = 0.2,
    cellsize = 1.0,
    circular = false,
    adjust = false,
    erode = false,
)
    _pmf2(A, ωₘ, slope, dhₘ, dh₀, cellsize, circular, adjust, erode)
end

function _pmf2(
    A::AbstractMatrix{<:Real},
    ωₘ::AbstractMatrix{<:Real},
    slope::Real,
    dhₘ::Real,
    dh₀::Real,
    cellsize::Real,
    circular::Bool,
    adjust::Bool,
    erode::Bool,
)
    _pmf(
        A,
        ωₘ,
        Fill(slope, size(A)),
        dhₘ,
        Fill(dh₀, size(A)),
        cellsize,
        circular,
        adjust,
        erode,
    )
end

"""
    round_odd(x)

Rounds `x` to the nearest odd number.
"""
round_odd(x) = round(Int, x / 2, RoundDown) * 2 + 1

function halve_range(x, stop = 3)
    out = [x]
    while x > stop
        x = round_odd(x / 2)
        insert!(out, 1, x)
    end
    return out
end

function _pmf2(
    A::AbstractMatrix{<:Real},
    windows::AbstractMatrix{<:Real},
    slope::AbstractMatrix{<:Real},
    dhₘ::Real,
    dh₀::AbstractMatrix{<:Real},
    cellsize::Real,
    circular::Bool,
    adjust::Bool,
    erode::Bool,
)

    # Compute windowsizes and thresholds
    iwindows = round.(Int, windows ./ cellsize) .+ 3
    ωₘ = maximum(iwindows)
    # @info ωₘ
    # κ_max = floor(Int, log2(ωₘ - 1))  # determine # iterations based on exp growth
    # windowsizes = Int.(exp2.(1:κ_max)) .+ 1
    windowsizes = halve_range(ωₘ)
    # @info windowsizes

    # Compute tresholds
    dwindows = vcat(windowsizes[1], windowsizes)  # prepend first element so we get 0 as diff
    window_diffs = [dwindows[i] - dwindows[i - 1] for i in 2:length(dwindows)]
    # height_tresholds = [min(dhₘ, slope * window_diff * cellsize + dh₀) for window_diff in window_diffs]

    # Set up arrays
    Af = copy(A)  # array to be morphed
    nan_mask = isnan.(Af)
    Af[nan_mask] .= Inf  # Replace NaN with Inf, as to always filter these

    B = copy(A)  # max_elevation raster
    out = copy(A)  # max_elevation raster

    flags = similar(A, Float64)  # 0 = ground, other values indicate window size
    fill!(flags, 0.0)
    flags[nan_mask] .= NaN

    mask = falses(size(A))

    # Iterate over window sizes and height tresholds
    for (i, ωₖ) in enumerate(windowsizes)
        # @info i, ωₖ
        s = (i > 1) && adjust ? dilate(slope, window_diffs[i]) : slope
        @debug "Window $(ωₖ), $(window_diffs[i]) slope sum: $(sum(s))"
        dhₜ = min.(dhₘ, s * window_diffs[i] * cellsize .+ dh₀)
        if erode
            if circular
                mapwindowcirc_approx!(minimum_mask, B, ωₖ, Af, Inf)
            else
                # mapwindow_stack!(minimum, A, ωₖ, Af)
                Af = LocalFilters.erode(out, ωₖ)
            end
        else
            if circular
                opening_circ!(Af, ωₖ, B)
            else
                Af = LocalFilters.opening(out, ωₖ)
            end
        end
        mask .= (out .- Af) .> dhₜ
        mask .&= (ωₖ .<= iwindows)
        for I in eachindex(flags)
            if mask[I]
                flags[I] = ωₖ
            end
        end
        out = Af
        # B[mask] .= min.(B[mask], Af[mask] .+ dhₜ[mask])
    end
    Af, flags
end
