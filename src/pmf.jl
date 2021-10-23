"""
```
B, flags = pmf(A; ωₘ, slope, dhₘ, dh₀, cellsize)
```
Applies the progressive morphological filter by *Zhang et al. (2003)* [^zhang2003] to `A`.

# Output
- `B::Array{T,2}` Maximum allowable values
- `flags::Array{Float64,2}` A sized array with window sizes if filtered, zero if not filtered.

Afterwards, one can retrieve the resulting mask for `A` by `A .<= B` or `flags .== 0.`.

# Arguments
- `A::Array{T,2}` Input Array
- `ωₘ::Float64=20.` Maximum window size [m]
- `slope::Float64=0.01` Terrain slope [m/m]
- `dhₘ::Float64=2.5` Maximum elevation threshold [m]
- `dh₀::Float64=0.2` Initial elevation threshold [m]
- `cellsize::Float64=1.` Cellsize in [m]

[^zhang2003]: Zhang, Keqi, Shu-Ching Chen, Dean Whitman, Mei-Ling Shyu, Jianhua Yan, and Chengcui Zhang. “A Progressive Morphological Filter for Removing Nonground Measurements from Airborne LIDAR Data.” IEEE Transactions on Geoscience and Remote Sensing 41, no. 4 (2003): 872–82. <https://doi.org/10.1109/TGRS.2003.810682>.
"""
function pmf(A::Array{T,2};
    ωₘ::Float64=20.,
    slope::Float64=0.01,
    dhₘ::Float64=2.5,
    dh₀::Float64=0.2,
    cellsize::Float64=1.0) where T <: Real

    # Compute windowsizes and thresholds
    ωₘ = round(Int, ωₘ / cellsize)
    κ_max = floor(Int, log2(ωₘ - 1))  # determine # iterations based on exp growth
    windowsizes = Int.(exp2.(1:κ_max)) .+ 1

    # Compute tresholds
    dwindows = vcat(windowsizes[1], windowsizes)  # prepend first element so we get 0 as diff
    window_diffs = [dwindows[i] - dwindows[i - 1] for i in 2:length(dwindows)]
    height_tresholds = [min(dhₘ, slope * window_diff * cellsize + dh₀) for window_diff in window_diffs]

    # Set up arrays
    Af = copy(A)  # array to be morphed
    nan_mask = isnan.(Af)
    Af[nan_mask] .= Inf  # Replace NaN with Inf, as to always filter these

    B = copy(A)  # max_elevation raster
    out = copy(A)  # max_elevation raster

    flags = zeros(size(A))  # 0 = ground, other values indicate window size
    flags[nan_mask] .= NaN

    mask = falses(size(A))

    # Iterate over window sizes and height tresholds
    for (ωₖ, dhₜ) in zip(windowsizes, height_tresholds)
        opening!(Af, ωₖ, out)
        mask .= (A .- Af) .> dhₜ
        for I in eachindex(flags)
            if mask[I] && flags[I] == 0
                flags[I] = ωₖ
            end
        end
        B .= min.(B, Af .+ dhₜ)
    end

    B, flags
end
