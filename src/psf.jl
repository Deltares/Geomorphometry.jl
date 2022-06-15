"""
```
B, flags = psf(A; ωₘ, slope, dhₘ, dh₀, cellsize)
```
Applies a progressive slope filter to `A`.

# Output
- `B::Array{T,2}` Maximum allowable values based A + slope, dhₘ, dh₀
- `flags::Array{Float64,2}` A sized array with window sizes if filtered, zero if not filtered.

Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.

# Arguments
- `A::Array{T,2}` Input Array
- `ωₘ::Float64=20.` Maximum window size [m]
- `slope::Float64=0.01` Terrain slope [m/m]
- `dhₘ::Float64=2.5` Maximum elevation threshold [m]
- `dh₀::Float64=0.2` Initial elevation threshold [m]
- `cellsize::Float64=1.` Cellsize in [m]
"""
function apsf(A::AbstractMatrix{T};
    ωₘ::Float64=20.0,
    slope::Matrix{Float64},
    dhₘ::Float64=2.5,
    dh₀::Float64=0.2,
    cellsize::Float64=1.0,
    circular=false) where {T<:Real}

    # Compute windowsizes and thresholds
    ωₘ = round(Int, ωₘ / cellsize)
    κ_max = floor(Int, log2(ωₘ - 1))  # determine # iterations based on exp growth
    windowsizes = Int.(exp2.(1:κ_max)) .+ 1

    # Compute tresholds
    dwindows = vcat(windowsizes[1], windowsizes)  # prepend first element so we get 0 as diff
    window_diffs = [dwindows[i] - dwindows[i-1] for i = 2:length(dwindows)]
    # height_tresholds = [min(dhₘ, slope[1] * window_diff * cellsize + dh₀) for window_diff in window_diffs]
    # @debug "Using the following thresholds: $height_tresholds for the following windows: $windowsizes"

    # Set up arrays
    Af = copy(A)  # array to be morphed
    nan_mask = isnan.(Af)
    Af[nan_mask] .= Inf  # Replace NaN with Inf, as to always filter these

    B = copy(A)  # max_elevation raster

    flags = zeros(size(A))  # 0 = ground, other values indicate window size
    flags[nan_mask] .= NaN

    mask = falses(size(A))
    dhₜ = min.(dhₘ, slope * window_diffs[1] * cellsize .+ dh₀)

    # Iterate over window sizes and height tresholds
    for (i, ωₖ) in enumerate(windowsizes)
        dhₜ = min.(dhₘ, slope * window_diffs[i] * cellsize .+ dh₀)
        if circular
            mapwindowcirc_approx!(minimum_mask, A, ωₖ, Af, Inf)
        else
            mapwindow_stack!(minimum, A, ωₖ, Af)
        end
        mask .= (A .- Af) .> dhₜ
        for I in eachindex(flags)
            if mask[I] && flags[I] == 0
                flags[I] = ωₖ
            end
        end
        B .= min.(B, Af .+ dhₜ)
    end

    B, flags, dhₜ
end



function psf(A::AbstractMatrix{T};
    ωₘ::Float64=20.0,
    slope::Float64=0.01,
    dhₘ::Float64=2.5,
    dh₀::Float64=0.2,
    cellsize::Float64=1.0,
    circular=false) where {T<:Real}
    return apsf(A;
        ωₘ=ωₘ,
        slope=fill(slope, size(A)),
        dhₘ=dhₘ,
        dh₀=dh₀,
        cellsize=cellsize,
        circular=circular)
end
