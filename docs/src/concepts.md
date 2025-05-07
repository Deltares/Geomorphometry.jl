# Concepts

In Geomorphometry.jl we provide a set of tools to analyse and visualize the shape of the Earth.
The package is designed to be fast, flexible, and easy to use. It is built with the following concepts in mind:

## Geospatially aware
With package extensions on [GeoArrays.jl](https://github.com/evetion/GeoArrays.jl) and [Rasters.jl](https://github.com/rafaqz/Rasters.jl) geospatial data is automatically handled to set the correct cellsize, even for geographical DEMs.

```@setup plots
using Geomorphometry, CairoMakie, GeoArrays, Rasters
set_theme!(theme_minimal(); transparency = true)
CairoMakie.activate!(type = "png")
r = GeoArrays.read("saba.tif")
mask = ismissing.(r)
dtm = coalesce(r, NaN)
```

:::tabs

== Rasters (projected)
```@example plots
r = Raster("saba.tif")
Geomorphometry.cellsize(r)
```
```@example plots
heatmap(multihillshade(r))
```

== GeoArrays (projected)
```@example plots
r = GeoArrays.read("saba.tif")
r = coalesce(r, NaN)
Geomorphometry.cellsize(r)
```
```@example plots
heatmap(multihillshade(r))
```

== Rasters (geographic)
```@example plots
r = Raster("Copernicus_DSM_10_N52_00_E004_00_DEM.tif")
Geomorphometry.cellsize(r)
```
```@example plots
heatmap(multihillshade(r))
```

== GeoArrays (geographic)
```@example plots
r = GeoArrays.read("Copernicus_DSM_10_N52_00_E004_00_DEM.tif")
r = coalesce(r, NaN)
Geomorphometry.cellsize(r)
```
```@example plots
heatmap(multihillshade(r))
```

:::

## Multiple algorithms
We have implemented several algorithms for a common operation so that you can choose the one that best fits your needs. For example, the `slope` function has three methods: `Horn`, `ZevenbergenThorne`, and `MaximumDownwardGradient`, as shown in the [Usage](usage.md) section.

Sometimes, as is the case for the `FD8` algorithm, these methods take different parameters that influence the results. `FD8` takes a `p` parameter that is used to weigh the flow direction, with higher powers resulting in less divergent flows (and thus more like D8).

:::tabs

== FD8 with power of 1
```@example plots
acc, ldd = flowaccumulation(dtm; method=FD8(1))
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== FD8 with power of 2
```@example plots
acc, ldd = flowaccumulation(dtm; method=FD8(2))
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== FD8 with power of 5
```@example plots
acc, ldd = flowaccumulation(dtm; method=FD8(5))
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```
== D8 (single flow direction)
```@example plots
acc, ldd = flowaccumulation(dtm; method=D8())
acc[mask] .= NaN  # hide
heatmap(log10.(acc); colormap=:rain)
```

:::


## Multiple scales
Inspired by the excellent [MultiScaleDTM](https://github.com/ailich/MultiscaleDTM) package in R by [ilichMultiscaleDTMOpensourcePackage2023](@citet), we have added multiscale options to some filters.

Relative terrain filters have a `window` keyword argument for a Stencil from [Stencils.jl](https://github.com/rafaqz/Stencils.jl) package.

:::tabs

== Square window of 1
```@example plots
Geomorphometry.Moore(1)
```
```@example plots
heatmap(TPI(dtm, Geomorphometry.Moore(1)); colorrange=(0,25), colormap=:speed)
```

== Square window of 3
```@example plots
Geomorphometry.Moore(3)
```
```@example plots
heatmap(TPI(dtm, Geomorphometry.Moore(3)); colorrange=(0,25), colormap=:speed)
```

== Square window of 5
```@example plots
Geomorphometry.Moore(5)
```
```@example plots
heatmap(TPI(dtm, Geomorphometry.Moore(5)); colorrange=(0,25), colormap=:speed)
```

:::


Other methods that require a specific type of window now take a `radius` kwarg, scaling said window.

:::tabs

== Radius of 1
```@example plots
Geomorphometry.scaled8nb(1)
```
```@example plots
heatmap(profile_curvature(dtm, radius=1); colorrange=(-1,1), colormap=:tarn)
```
== Radius of 3
```@example plots
Geomorphometry.scaled8nb(3)
```
```@example plots
heatmap(profile_curvature(dtm, radius=3); colorrange=(-1,1), colormap=:tarn)
```
== Radius of 5
```@example plots
Geomorphometry.scaled8nb(5)
```
```@example plots
heatmap(profile_curvature(dtm, radius=5); colorrange=(-1,1), colormap=:tarn)
```

:::
