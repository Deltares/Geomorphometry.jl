
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

![](pzzoijk.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(hillshade(ndem))
```

![](wankxuv.png){width=600px height=450px}

== GDAL

```julia
heatmap(hillshade(Geomorphometry.GDAL(), dem))
```

![](olpbpif.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(hillshade(Geomorphometry.GDAL(), ndem))
```

![](glggomt.png){width=600px height=450px}

:::

## Multihillshade {#Multihillshade}

:::tabs

== Ours

```julia
heatmap(multihillshade(dem))
```

![](frxwgav.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(multihillshade(ndem))
```

![](dwtffbz.png){width=600px height=450px}

== GDAL

```julia
heatmap(multihillshade(Geomorphometry.GDAL(), dem))
```

![](jpxyple.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(multihillshade(Geomorphometry.GDAL(), ndem))
```

![](ohevpje.png){width=600px height=450px}

:::

## Slope {#Slope}

### Horn {#Horn}

:::tabs

== Ours

```julia
heatmap(slope(dem, method=Horn()))
```

![](vejnovs.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(slope(ndem, method=Horn()))
```

![](gmhlfsz.png){width=600px height=450px}

== GDAL

```julia
heatmap(slope(Geomorphometry.GDAL(), dem, method=Horn()))
```

![](djpsgur.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(slope(Geomorphometry.GDAL(), ndem, method=Horn()))
```

![](ngogwqp.png){width=600px height=450px}

:::

### ZevenbergenThorne {#ZevenbergenThorne}

:::tabs

== Ours

```julia
heatmap(slope(dem, method=ZevenbergenThorne()))
```

![](ewrwmln.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(slope(ndem, method=ZevenbergenThorne()))
```

![](oxvewkk.png){width=600px height=450px}

== GDAL

```julia
heatmap(slope(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```

![](xlyjcrr.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(slope(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

![](hzkwkkc.png){width=600px height=450px}

:::

## Aspect {#Aspect}

### Horn {#Horn-2}

:::tabs

== Ours

```julia
heatmap(aspect(dem, method=Horn()))
```

![](xlbegkt.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(aspect(ndem, method=Horn()))
```

![](ihxwxjy.png){width=600px height=450px}

== GDAL

```julia
heatmap(aspect(Geomorphometry.GDAL(), dem, method=Horn()))
```

![](xopkxpj.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=Horn()))
```

![](oypnecj.png){width=600px height=450px}

:::

### ZevenbergenThorne {#ZevenbergenThorne-2}

:::tabs

== Ours

```julia
heatmap(aspect(dem, method=ZevenbergenThorne()))
```

![](pjydxjx.png){width=600px height=450px}

== Ours (reverse y-axis)

```julia
heatmap(aspect(ndem, method=ZevenbergenThorne()))
```

![](gsjwmuo.png){width=600px height=450px}

== GDAL

```julia
heatmap(aspect(Geomorphometry.GDAL(), dem, method=ZevenbergenThorne()))
```

![](humlfgp.png){width=600px height=450px}

== GDAL (reverse y-axis)

```julia
heatmap(aspect(Geomorphometry.GDAL(), ndem, method=ZevenbergenThorne()))
```

![](hcpdchm.png){width=600px height=450px}

:::
