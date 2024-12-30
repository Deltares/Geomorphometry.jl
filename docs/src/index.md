```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Geomorphometry.jl"
  tagline: "Analyzing and visualizing the shape of the Earth"
  image:
    src: logo.svg
    alt: Geomorphometry
  actions:
    - theme: brand
      text: Get Started
      link: /tutorials/installation.md
    - theme: alt
      text: View on Github
      link: https://github.com/Deltares/Geomorphometry.jl
    - theme: alt
      text: API Reference
      link: /reference
      
features:
  - title: Common operations
    details: Defines common methods for reading, writing and manipulating geospatial vector data, as one-liners, and as fast as possible.
    link: /tutorials/usage
  - title: Multiple algorithms
    details: Supports hillshade, multihillshade and perceptual shading for visualisation of elevation models.
    link: /background/formats
  - title: Seamless integration
    details: Geomorphometry.jl is fully compatible with the AbstractArray and GeoInterface.jl ecosystems. This enables plotting, operations and analysis using the full power of the Julia ecosystem.


---
```

```@meta
CurrentModule = Geomorphometry
```

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
