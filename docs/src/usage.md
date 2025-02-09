# Usage

```@meta
CurrentModule = Geomorphometry
```

```@setup plots
using Geomorphometry, CairoMakie, GeoArrays, Statistics # hide
set_theme!(theme_minimal(); transparency = true) # hide
CairoMakie.activate!(type = "png")
A = GeoArrays.read("saba.tif") # hide
mask = ismissing.(A) # hide
dem = coalesce(A, NaN).A # hide
```

In Geomorphometry.jl we provide a set of tools to analyze and visualize the shape of the Earth. The package is designed to be fast, flexible, and easy to use.
Moreover, we have implemented several algorithms for a common operation so that you can choose the one that best fits your needs.

In these pages we will use the elevation model of [Saba](https://en.wikipedia.org/wiki/Saba_(island)) to showcase the different categories of operations that are available in Geomorphometry.jl.

## Visualization
Visualization is done using the [`hillshade`](@ref), [`multihillshade`](@ref), and [`pssm`](@ref) functions. The first two shade the terrain by using a single or multiple light source(s) respectively, while `pssm` is a slope map exaggerated for human perception.

```@example plots
s = multihillshade(dem; cellsize=5)
# s[mask] .= NaN  # hide
heatmap(s)
```

We can also overlay the `pssm` visualization on top of a colored `aspect` map.

```@example plots
f = heatmap(aspect(dem); colormap=:curl)
heatmap!(pssm(dem, cellsize=5); colormap=Reverse(:greys), alpha=0.5)
f
```

## Derivatives
Common derivatives are implemented in Geomorphometry.jl. These include [`slope`](@ref), [`aspect`](@ref), and `curvature`. The latter is ill-defined, here we provide [`plan_curvature`](@ref) (also called *projected contour curvature*), [`profile_curvature`](@ref) (also called *normal slope line curvature*), and [`tangential_curvature`](@ref) (also called *normal contour curvature*). Note that functions here allow for a custom radius (but fixed positions, see X), as demonstrated for `profile_curvature`.

:::tabs

== Slope
```@example plots
heatmap(slope(dem, cellsize=5); colormap=:matter)
```
== Aspect
```@example plots
heatmap(aspect(dem); colormap=:phase)
```
== Curvature
```@example plots
heatmap(profile_curvature(dem, cellsize=5, radius=3); colorrange=(-1,1), colormap=:tarn)
```

:::

## Relative position
There are several terrain descriptors that can be used to analyze the relative position of a point with respect to its neighbors. These include [`TPI`](@ref), [`TRI`](@ref), [`RIE`](@ref), [`BPI`](@ref),and [`roughness`](@ref). Here we use `BPI`, but with a custom sized `Window` (from Stencils.jl). All these functions can be used with a custom window size.

:::tabs

== BPI
```@example plots
heatmap(BPI(dem, Geomorphometry.Annulus(5, 3)); colormap=:delta, colorrange=(-10,10))
```

:::

## Hydrology
Hydrological operations are used to analyze the flow of water on the terrain. We provide [`filldepressions`](@ref) to fill depressions, and [`flowaccumulation`](@ref) to calculate the flow accumulation. Here we use `flowaccumulation` to calculate the flow accumulation. Note that the local drainage direction is also returned. By default the FD8 algorithm is used, but you can also use the D∞ or D8 algorithm by setting the `method` keyword argument to `DInf()` or `D8()`.

:::tabs

== Flow accumulation with FD8
```@example plots
acc, ldd = flowaccumulation(dem; method=FD8(), cellsize=5)
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D∞
```@example plots
acc, ldd = flowaccumulation(dem; method=DInf(), cellsize=5)
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Flow accumulation with D8
```@example plots
acc, ldd = flowaccumulation(dem; method=D8(), cellsize=5)
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== Underlying method

We use the Priority Flood method by Barnes et al. (2014).

```@example plots
out = zeros(size(dem))
queued = falses(size(dem))
fig, ax, hm = heatmap(out, colormap=:turbo, colorrange=(0, 600))

openn = Geomorphometry.PriorityQueue{CartesianIndex{2}, eltype(dem)}()
pit = Geomorphometry.DataStructures.Queue{CartesianIndex{2}}()

R = CartesianIndices(dem)
I_first, I_last = first(R), last(R)
fps=60
for cell in Geomorphometry.edges(R)
    Geomorphometry.enqueue!(openn, cell, dem[cell])
    queued[cell] = true  # queued
end
record(fig, expanduser("test.mp4")) do io
    i = 0
    while (!isempty(openn) || !isempty(pit))
        cell = !isempty(pit) ? Geomorphometry.DataStructures.dequeue!(pit) : Geomorphometry.dequeue!(openn)
        out[cell] = dem[cell]
        # sleep(1/fps)
        # hm[1] = out # update data
        if iszero(i % 5_000)
            hm[1] = out # update data
            recordframe!(io)
        end
        for ncell in max(I_first, cell - Geomorphometry.Δ):min(I_last, cell + Geomorphometry.Δ)
            (queued[ncell] || ncell == cell) && continue
            queued[ncell] = true
            if dem[ncell] <= dem[cell]
                dem[ncell] = dem[cell]
                Geomorphometry.DataStructures.enqueue!(pit, ncell)
            else
                Geomorphometry.enqueue!(openn, ncell, dem[ncell])
            end
        end
        i += 1
    end
end
```
```@raw html
<video autoplay loop muted playsinline controls src="./test.mp4" />
```

:::

Combined with the previous analysis, we can calculate [`TWI`](@ref) and [`SPI`](@ref) indices. These indices are used to analyze the terrain's ability to accumulate water and the terrain's ability to store water respectively.

```@example plots
twi = TWI(dem; method=D8(), cellsize=5)
heatmap(twi; colormap=:tempo)
```

```@example plots
# labels = filldepressions(dem[250:750, 250:750])[2]
labels = filldepressions(dem)[2]
# GeoArrays.write("labels.tif", GeoArray(Int32.(labels)))  # hide
random_colormap = [Makie.ColorTypes.RGB(rand(), rand(), rand()) for _ in 1:length(unique(labels))]
# heatmap(labels; colormap=cgrad(:rainbow))
```


## Terrain filters
While filters seem unrelated to the previously discussed methods, these are used to filter (rasterized) pointclouds i.e. surface models (DSM) to arrive at a terrain model (DTM). We provide [`pmf`](@ref), [`smf`](@ref), and [`skb`](@ref) filters. Here we use the `pmf` filter to remove elevations that do not pass the slope filter (normally vegetation and buildings).

```@example plots
B, flags = pmf(dem, ωₘ = 20.0, slope = 0.11, dhₘ = 2.5, dh₀ = 0.2, cellsize = 5.0)
A = copy(dem)
A[A .> B] .= NaN
heatmap(A)
```


## Cost distance
Finally there are so-called cost distance functions. Instead of calculating a distance in the Euclidean sense, these functions calculate the cost of moving from one point to another, given a spatial varying cost. We provide a single [`spread`](@ref) function with multiple `method`s. The function takes starting points, an initial cost, and a cost (friction) map.

```@example plots
# What's the shortest path to the coast from the top, taking into account the steepness of your path?
dem[isnan.(dem)].=500  # hide
dist, zone = spread(dem .> 850, 0, slope(dem, cellsize=5).^2, method=FastSweeping())
h = contour(copy(dist))
dist[dem .> 1] .= Inf  # only consider coastal cells
heatmap!(pssm(dem, cellsize=5); colormap=Reverse(:greys), alpha=0.25)
scatter!(Tuple(findmin(dist)[2]), color=:red)
h
```
