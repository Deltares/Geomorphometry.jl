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
    details: Defines common methods for analyzing, filtering and visualizing (global) elevation models. All methods are implemented in Julia and are fast and scalable.
    link: /usage
  - title: Multiple algorithms
    details: Choose from multiple algorithms, ranging from different derivations of slope to multiresolution stencils. This enables the user to select the most appropriate method for their use case.
    link: /concepts
  - title: Seamless integration
    details: Geomorphometry.jl is fully compatible with the AbstractArray and GeoInterface.jl ecosystems. This enables plotting, operations and analysis using the full power of the Julia ecosystem.


---
```

```@meta
CurrentModule = Geomorphometry
```

A Julia package for researchers analysing the shape of the surface of the Earth ([Geomorphometry](https://en.wikipedia.org/wiki/Geomorphometry)) and its applications in hydrology, geomorphology, and environmental science.

## Functionality
The package provides efficient implementations of common geomorphometric operations, including:
- Terrain filters, such as Progressive Morphological Filters (PMF, SMF) and Skewness balancing
- Geospatial cost (friction) operations that mimic PCRaster. These functions should however be more Julian, extensible and scale better.
- Visualization, such as Perceptually Shaded Slope Map (PSSM)
- Terrain analysis functions, such as slope, aspect, roughness, Topographic Position Index (TPI), Terrain Ruggedness Index (TRI), curvature and hillslope.

When possible, multiple algorithms are provided for the same operation, allowing users to select the most appropriate method for their specific use case.

## Installation
The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Geomorphometry
```

## Publications
The code (specifically the slope filters) in this package were used to produce DeltaDTM, a global coastal digital terrain model (DTM) at 30 m resolution:
> Pronk, M., Hooijer, A., Eilander, D., Haag, A., de Jong, T., Vousdoukas, M., Vernimmen, R., Ledoux, H., & Eleveld, M. (2024). DeltaDTM: A global coastal digital terrain model. Scientific Data, 11(1), 273. https://doi.org/10.1038/s41597-024-03091-9

## Contributing
Contributions are welcome! Please feel free to report issues, feature requests, or submit pull requests on [GitHub](https://github.com/evetion/GeoDataFrames.jl/issues).
