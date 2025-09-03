
# Validation {#Validation}

This chapter compares our output with that of [gdaldem](https://gdal.org/en/stable/programs/gdaldem.html).

We test with the same elevation model as that in the [Usage](usage.md) page, and also include an upside down version of it, to ensure that the algorithms are robust to different axes conventions.

```julia
@info Geomorphometry.cellsize(dem)  # default negative y spacing
@info Geomorphometry.cellsize(ndem)  # reversed y spacing
```


```ansi
[91m[1m┌ [22m[39m[91m[1mError: [22m[39mUnknown CRS type, please report this issue for the given crs/file
[91m[1m└ [22m[39m[90m@ GeoArrays ~/.julia/packages/GeoArrays/X63Mm/src/geointerface.jl:11[39m
[36m[1m[ [22m[39m[36m[1mInfo: [22m[39m(5.0, -5.0)
[91m[1m┌ [22m[39m[91m[1mError: [22m[39mUnknown CRS type, please report this issue for the given crs/file
[91m[1m└ [22m[39m[90m@ GeoArrays ~/.julia/packages/GeoArrays/X63Mm/src/geointerface.jl:11[39m
[36m[1m[ [22m[39m[36m[1mInfo: [22m[39m(5.0, 5.0)
```


## Hillshade {#Hillshade}

:::tabs

== Ours

```julia
heatmap(hillshade(dem))
```

![](vvlqpkc.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(hillshade(ndem))
```

![](elusomt.png){width=600px height=450px}

== GDAL

```julia
heatmap(hillshade(Geomorphometry.GDAL(), dem))
```

![](hehjsqz.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(hillshade(Geomorphometry.GDAL(), ndem))
```

![](kyiwhar.png){width=600px height=450px}

:::

## Multihillshade {#Multihillshade}

:::tabs

== Ours

```julia
heatmap(multihillshade(dem))
```

![](evjtuoj.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(multihillshade(ndem))
```

![](goanwfz.png){width=600px height=450px}

== GDAL

```julia
heatmap(multihillshade(Geomorphometry.GDAL(), dem))
```

![](qsusfrt.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(multihillshade(Geomorphometry.GDAL(), ndem))
```

![](udirtfy.png){width=600px height=450px}

:::

## Slope {#Slope}

### Horn {#Horn}

:::tabs

== Ours

```julia
heatmap(slope(dem, method=Horn()))
```

![](hgyzczt.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(slope(ndem, method=Horn()))
```

![](zvpbeyu.png){width=600px height=450px}

== GDAL

```julia
heatmap(slope(Geomorphometry.GDAL(), dem, method=Horn()))
```

![](jwemclt.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(slope(Geomorphometry.GDAL(), ndem, method=Horn()))
```

![](gmkzigd.png){width=600px height=450px}

:::

### ZevenbergenThorne {#ZevenbergenThorne}

:::tabs

== Ours

```julia
heatmap(slope(dem, method=ZevenbergenThorne()))
```

![](rzugxfr.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(slope(ndem, method=ZevenbergenThorne()))
```

![](twarhcn.png){width=600px height=450px}

== GDAL

```julia
heatmap(slope(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```

![](bzynolz.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(slope(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

![](fvsywwm.png){width=600px height=450px}

:::

## Aspect {#Aspect}

### Horn {#Horn-2}

:::tabs

== Ours

```julia
heatmap(aspect(dem, method=Horn()))
```

![](jllmxmg.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(aspect(ndem, method=Horn()))
```

![](fyluiyf.png){width=600px height=450px}

== GDAL

```julia
heatmap(aspect(Geomorphometry.GDAL(), dem, method=Horn()))
```

![](uyqjznj.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=Horn()))
```

![](scwccuf.png){width=600px height=450px}

:::

### ZevenbergenThorne {#ZevenbergenThorne-2}

:::tabs

== Ours

```julia
heatmap(aspect(dem, method=ZevenbergenThorne()))
```

![](kvrtowi.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(aspect(ndem, method=ZevenbergenThorne()))
```

![](eiuvifj.png){width=600px height=450px}

== GDAL

```julia
heatmap(aspect(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```

![](ojdrskr.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

![](fxauraa.png){width=600px height=450px}

:::
