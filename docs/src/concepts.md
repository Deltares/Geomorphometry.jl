# Concepts

In Geomorphology.jl we provide a set of tools to analyse and visualize the shape of the Earth.
The package is designed to be fast, flexible, and easy to use.
Moreover, we have implemented several algorithms for a common operation so that you can choose the one that best fits your needs.

```@example
sleep(1.5)
```

```@example plots
using Geomorphometry, CairoMakie
A = rand(10,10)
s = slope(A)
heatmap(s)
```
