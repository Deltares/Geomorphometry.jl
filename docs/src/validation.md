# Validation

```@setup plots
using Geomorphometry, CairoMakie, GeoArrays, Rasters
set_theme!(theme_minimal(); transparency = true) 
CairoMakie.activate!(type = "png")
fn = "saba.tif"
A = GeoArrays.read(fn)
demr = Raster(fn; checkmem=false) 
dem = coalesce(A, NaN) 
ndem = GeoArrays.flipud!(deepcopy(dem))
# ndem = reverse(dem, dims=2)
```

This chapter compares our output with that of [gdaldem](https://gdal.org/en/stable/programs/gdaldem.html).

We test with the same elevation model as that in the [Usage](usage.md) page, and also include an upside down version of it, to ensure that the algorithms are robust to different axes conventions.

```@example plots
@info Geomorphometry.cellsize(dem)  # default negative y spacing
@info Geomorphometry.cellsize(ndem)  # reversed y spacing
```

## Hillshade

:::tabs

== Ours
```@example plots
heatmap(hillshade(dem))
```
== Ours (reverse y-axis)
```@example plots
heatmap(hillshade(ndem))
```
== GDAL
```@example plots
heatmap(hillshade(Geomorphometry.GDAL(), dem))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(hillshade(Geomorphometry.GDAL(), ndem))
```

:::

## Multihillshade

:::tabs

== Ours
```@example plots
heatmap(multihillshade(dem))
```
== Ours (reverse y-axis)
```@example plots
heatmap(multihillshade(ndem))
```
== GDAL
```@example plots
heatmap(multihillshade(Geomorphometry.GDAL(), dem))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(multihillshade(Geomorphometry.GDAL(), ndem))
```
:::

## Slope

### Horn
:::tabs

== Ours
```@example plots
heatmap(slope(dem, method=Horn()))
```
== Ours (reverse y-axis)
```@example plots
heatmap(slope(ndem, method=Horn()))
```
== GDAL
```@example plots
heatmap(slope(Geomorphometry.GDAL(), dem, method=Horn()))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(slope(Geomorphometry.GDAL(), ndem, method=Horn()))
```

:::

### ZevenbergenThorne

:::tabs

== Ours
```@example plots
heatmap(slope(dem, method=ZevenbergenThorne()))
```
== Ours (reverse y-axis)
```@example plots
heatmap(slope(ndem, method=ZevenbergenThorne()))
```
== GDAL
```@example plots
heatmap(slope(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(slope(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

:::


## Aspect

### Horn

:::tabs

== Ours
```@example plots
heatmap(aspect(dem, method=Horn()))
```
== Ours (reverse y-axis)
```@example plots
heatmap(aspect(ndem, method=Horn()))
```
== GDAL
```@example plots
heatmap(aspect(Geomorphometry.GDAL(), dem, method=Horn()))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=Horn()))
```

:::


### ZevenbergenThorne

:::tabs

== Ours
```@example plots
heatmap(aspect(dem, method=ZevenbergenThorne()))
```
== Ours (reverse y-axis)
```@example plots
heatmap(aspect(ndem, method=ZevenbergenThorne()))
```
== GDAL
```@example plots
heatmap(aspect(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```
== GDAL (reverse y-axis)
```@example plots
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

:::
