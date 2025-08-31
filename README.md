[![CI](https://github.com/Deltares/GeoRasterFiltering.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Deltares/GeoRasterFiltering.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/Deltares/Geomorphometry.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Deltares/Geomorphometry.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deltares.github.io/Geomorphometry.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltares.github.io/Geomorphometry.jl/dev/)

# Geomorphometry
A Julia package for researchers analysing the shape of the surface of the Earth ([Geomorphometry](https://en.wikipedia.org/wiki/Geomorphometry)) and its applications in hydrology, geomorphology, and environmental science.

## Functionality
The package provides efficient implementations of common geomorphometric operations, including:
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

## Alternative Packages
If you are working in Python the [xDEM](https://xdem.readthedocs.io/en/stable/) package provides a comprehensive suite of tools for DEM analysis
