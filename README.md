[![CI](https://github.com/Deltares/GeoRasterFiltering.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Deltares/GeoRasterFiltering.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/Deltares/Geomorphometry.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Deltares/Geomorphometry.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deltares.github.io/Geomorphometry.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltares.github.io/Geomorphometry.jl/dev/)

# Geomorphometry
Geospatial operations, cost and filtering algorithms as used in for elevation rasters.

## Functionality
- Terrain filters, such as Progressive Morphological Filters (PMF, SMF) and Skewness balancing
- Geospatial cost (friction) operations that mimic PCRaster. These functions should however be more Julian, extensible and scale better.
- Visualization, such as Perceptually Shaded Slope Map (PSSM)
- Terrain analysis functions, such as slope, aspect, roughness, Topographic Position Index (TPI), Terrain Ruggedness Index (TRI), curvature and hillslope.

## Installation
The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Geomorphometry
```
