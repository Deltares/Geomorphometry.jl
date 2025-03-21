# Usage

```@meta
CurrentModule = Geomorphometry
```

```@setup plots
using Geomorphometry, CairoMakie, GeoArrays, Statistics # hide
import Eikonal # hide
set_theme!(theme_minimal(); transparency = true) # hide
CairoMakie.activate!(type = "png")
A = GeoArrays.read("saba.tif") # hide
B = GeoArrays.read("saba_dsm.tif") # hide
mask = ismissing.(A) # hide
# dem = coalesce(A, NaN).A # hide
dtm = coalesce(A, NaN) # hide
dsm = coalesce(B, NaN) # hide
# ndem = GeoArrays.flipud!(deepcopy(dtm))

fdtm = deepcopy(dtm)
out = zeros(size(fdtm))
queued = falses(size(fdtm))
fig, ax, hm = heatmap(out, colormap=:turbo, colorrange=(0, 600))

openn = Geomorphometry.PriorityQueue{CartesianIndex{2}, eltype(fdtm)}()
pit = Geomorphometry.DataStructures.Queue{CartesianIndex{2}}()

R = CartesianIndices(fdtm)
I_first, I_last = first(R), last(R)
fps=60
for cell in Geomorphometry.edges(R)
    Geomorphometry.enqueue!(openn, cell, fdtm[cell])
    queued[cell] = true  # queued
end
record(fig, expanduser("test.mp4")) do io
    i = 0
    while (!isempty(openn) || !isempty(pit))
        cell = !isempty(pit) ? Geomorphometry.DataStructures.dequeue!(pit) : Geomorphometry.dequeue!(openn)
        out[cell] = fdtm[cell]
        # sleep(1/fps)
        # hm[1] = out # update data
        if iszero(i % 5_000)
            hm[1] = out # update data
            recordframe!(io)
        end
        for ncell in max(I_first, cell - Geomorphometry.Δ):min(I_last, cell + Geomorphometry.Δ)
            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true
            if fdtm[ncell] <= fdtm[cell]
                fdtm[ncell] = fdtm[cell]
                Geomorphometry.DataStructures.enqueue!(pit, ncell)
            else
                Geomorphometry.enqueue!(openn, ncell, fdtm[ncell])
            end
        end
        i += 1
    end
end
```

In Geomorphometry.jl we provide a set of tools to analyze and visualize the shape of the Earth. The package is designed to be fast, flexible, and easy to use.
Moreover, we have implemented several algorithms for a common operation so that you can choose the one that best fits your needs.

In these pages we will use the elevation model of [Saba](https://en.wikipedia.org/wiki/Saba_(island)) to showcase the different categories of operations that are available in Geomorphometry.jl.

## Visualization
Visualization is done using the [`hillshade`](@ref), [`multihillshade`](@ref), and [`pssm`](@ref) functions. The first two shade the terrain by using a single or multiple light source(s) respectively, while `pssm` is a slope map exaggerated for human perception.


:::tabs

== Hillshade
```@example plots
heatmap(hillshade(dtm))
```

== Multihillshade
```@example plots
heatmap(multihillshade(dtm, azimuth=0:60:270, zenith=60))
```

== PSSM
```@example plots
heatmap(pssm(dtm), colormap=Reverse(:viridis))
```
:::

We can also overlay the `pssm` visualization on top of a colored `aspect` map.

```@example plots
f = heatmap(aspect(dtm); colormap=:curl)
heatmap!(pssm(dtm); colormap=Reverse(:greys), alpha=0.5)
f
```

## Derivatives
Common derivatives are implemented in Geomorphometry.jl. These include [`slope`](@ref), [`aspect`](@ref), and `curvature`. The latter is ill-defined, here we provide [`plan_curvature`](@ref) (also called *projected contour curvature*), [`profile_curvature`](@ref) (also called *normal slope line curvature*), and [`tangential_curvature`](@ref) (also called *normal contour curvature*). Note that functions here allow for a custom radius (but fixed positions, see X), as demonstrated for `profile_curvature`.


### Slope

:::tabs

== Horn
```@example plots
heatmap(slope(dtm; method=Horn()); colormap=:matter, colorrange=(0, 60))
```
== ZevenbergenThorne
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne()); colormap=:matter, colorrange=(0, 60))
```
== MaximumDownwardGradient
```@example plots
heatmap(slope(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:matter, colorrange=(0, 60))
```

:::

We can also determine the slope of the terrain along a specific direction.

:::tabs

== Horn (0°)
```@example plots
heatmap(slope(dtm; method=Horn(), direction=0); colormap=:matter, colorrange=(-45, 45))
```
== ZevenbergenThorne (0°)
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=0); colormap=:matter, colorrange=(-45, 45))
```
== Horn (90°)
```@example plots
heatmap(slope(dtm; method=Horn(), direction=90); colormap=:matter, colorrange=(-45, 45))
```
== ZevenbergenThorne (90°))
```@example plots
heatmap(slope(dtm; method=ZevenbergenThorne(), direction=90); colormap=:matter, colorrange=(-45, 45))
```

:::

### Aspect

:::tabs

== Horn
```@example plots
heatmap(aspect(dtm; method=Horn()); colormap=:romaO)
```
== ZevenbergenThorne
```@example plots
heatmap(aspect(dtm; method=ZevenbergenThorne()); colormap=:romaO)
```
== MaximumDownwardGradient
```@example plots
heatmap(aspect(dtm; method=Geomorphometry.MaximumDownwardGradient()); colormap=:romaO)
```

:::

### Curvature

:::tabs

== Profile curvature
```@example plots
heatmap(profile_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

== Plan curvature
```@example plots
heatmap(plan_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

== Tangential curvature
```@example plots
heatmap(tangential_curvature(dtm); colorrange=(-1,1), colormap=:tarn)
```

:::

We can also determine the curvature of the terrain along a specific direction.

:::tabs

== Profile curvature in 90° direction
```@example plots
heatmap(profile_curvature(dtm, radius=1, direction=90); colorrange=(-1,1), colormap=:tarn)
```

== Plan curvature in 90° direction
```@example plots
heatmap(plan_curvature(dtm, radius=1, direction=90); colorrange=(-1,1), colormap=:tarn)
```

== Tangential curvature in 90° direction
```@example plots
heatmap(tangential_curvature(dtm, radius=1, direction=90); colorrange=(-1,1), colormap=:tarn)
```

== Profile curvature in 0° direction
```@example plots
heatmap(profile_curvature(dtm, radius=1, direction=0); colorrange=(-1,1), colormap=:tarn)
```

== Plan curvature in 0° direction
```@example plots
heatmap(plan_curvature(dtm, radius=1, direction=0); colorrange=(-1,1), colormap=:tarn)
```

== Tangential curvature in 0° direction
```@example plots
heatmap(tangential_curvature(dtm, radius=1, direction=0); colorrange=(-1,1), colormap=:tarn)
```

:::

## Relative position
There are several terrain descriptors that can be used to analyze the relative position of a point with respect to its neighbors. These include [`TPI`](@ref), [`TRI`](@ref), [`RIE`](@ref), [`BPI`](@ref), [`rugosity`](@ref) and [`roughness`](@ref). Here we use `BPI`, but with a custom sized `Window` (from Stencils.jl). All these functions can be used with a custom window size.

:::tabs

== BPI
```@example plots
heatmap(BPI(dtm, Geomorphometry.Annulus(5, 3)); colormap=:delta, colorrange=(-10,10))
```
== TPI
```@example plots
heatmap(TPI(dtm); colormap=:delta, colorrange=(-10,10))
```
== TRI
```@example plots
heatmap(TRI(dtm); colormap=:speed)
```
== RIE
```@example plots
heatmap(RIE(dtm); colormap=:speed)
```
== rugosity
```@example plots
heatmap(rugosity(dtm); colormap=:speed)
```
== roughness
```@example plots
heatmap(roughness(dtm); colormap=:speed)
```

:::

## Hydrology
Hydrological operations are used to analyze the flow of water on the terrain. We provide [`filldepressions`](@ref) to fill depressions, and [`flowaccumulation`](@ref) to calculate the flow accumulation. Here we use `flowaccumulation` to calculate the flow accumulation. Note that the local drainage direction is also returned. By default the FD8 algorithm is used, but you can also use the D∞ or D8 algorithm by setting the `method` keyword argument to `DInf()` or `D8()`.

:::tabs

== Flow accumulation with FD8
```@example plots
acc, ldd = flowaccumulation(dtm; method=FD8(2))
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D∞
```@example plots
acc, ldd = flowaccumulation(dtm; method=DInf())
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D8
```@example plots
acc, ldd = flowaccumulation(dtm; method=D8())
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Underlying method

We use the Priority Flood method by [barnesPriorityFloodOptimalDepressionFilling2014](@citet).

```@raw html
<video autoplay loop muted playsinline controls src="./test.mp4" />
```

:::

Combined with the previous analysis, we can calculate [`TWI`](@ref) and [`SPI`](@ref) indices. These indices are used to analyze the terrain's ability to accumulate water and the terrain's ability to store water respectively.


:::tabs

== TWI
```@example plots
twi = TWI(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```
== SPI
```@example plots
twi = SPI(dtm; method=FD8())
heatmap(twi; colormap=:tempo)
```

:::


## Terrain filters
While filters seem unrelated to the previously discussed methods, these are used to filter (rasterized) pointclouds i.e. surface models (DSM) to arrive at a terrain model (DTM). We provide [`pmf`](@ref), [`smf`](@ref), and [`skb`](@ref) filters. Here we use the `pmf` filter to remove elevations that do not pass the slope filter (normally vegetation and buildings).


:::tabs

== PMF
```@example plots
B, flags = pmf(dsm, ωₘ = 15.0, slope = 0.2, dhₘ = 5.0, dh₀ = 0.3)
A = copy(dsm)
A[A .> B] .= NaN
heatmap(A)
```
== SMF
```@example plots
A = smf(dsm; ω=15.0, slope=0.2)
heatmap(A)
```
== SKB
```@example plots
B, flags = skb(dsm)
A = copy(dsm)
A[A .> B] .= NaN
heatmap(A)
```
:::


## Cost distance
Finally there are so-called cost distance functions. Instead of calculating a distance in the Euclidean sense, these functions calculate the cost of moving from one point to another, given a spatial varying cost. We provide a single [`spread`](@ref) function with multiple `method`s. The function takes starting points, an initial cost, and a cost (friction) map.

Here we use the top of the volcano as the starting point, and the squared slope as the cost map. We compare the `Tomlin` method with the `Eikonal` and `Eastman` methods. `Eastman` requires a high number of iterations to converge (even at 25 it's not enough), and is thus slower.

:::tabs

== Tomlin
```@example plots
heatmap(spread(dtm .> 850, 0, slope(dtm).^2), colormap=:dense)
```
== Eikonal
```@example plots
heatmap(spread(dtm .> 850, 0, slope(dtm).^2, method=FastSweeping()), colormap=:dense)
```
== Eastman
```@example plots
heatmap(spread(dtm .> 850, 0, slope(dtm).^2, method=Eastman(iterations=25)), colormap=:dense)
```
:::
