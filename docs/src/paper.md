---
title: 'Geomorphometry.jl: Analyzing and visualizing the shape of the Earth in Julia.'
tags:
  - Julia
  - geomorphometry
  - GIS
  - hydrology
  - digital elevation model
  - terrain analysis
  - LIDAR
authors:
  - name: Maarten Pronk
    orcid: 0000-0001-8758-3939
    corresponding: true
    affiliation: 1
  - name: Rafael Schouten
    orcid: 0000-0002-8380-0884
    affiliation: 2
affiliations:
 - name: Deltares, Delft, The Netherlands
   index: 1
   ror: 01deh9c76
 - name: University of Melbourne, Melbourne, Australia
   index: 2
   ror: 01ej9dk98
date: 7 March 2026
bibliography: refs.bib

---

# Summary

Digital elevation models (DEMs) are fundamental datasets in the earth sciences, used in applications from flood risk assessment to habitat mapping and climate modelling.
Analysing the shape of the Earth's surface---the domain of geomorphometry---requires algorithms for deriving terrain attributes (slope, aspect, curvature), classifying landforms, routing water flows, and visualising terrain.
`Geomorphometry.jl` is a Julia [@bezansonJuliaFreshApproach2017] package that provides efficient, composable, and extendable implementations of these algorithms, designed to work with plain matrices as well as geospatial data types from the Julia ecosystem.

# Statement of need

Of the three widely used dynamic languages---Julia, Python, and R---Julia was the only one lacking a dedicated package for geomorphometric analysis.
`Geomorphometry.jl` fills this gap, enabling users to combine terrain analysis with Julia's broader scientific computing stack---differential equations, optimisation, GPU computing, and machine learning---without leaving the language.

The package implements over 30 functions spanning several domains (non-exhaustive):

- **Terrain derivatives**: slope, aspect, and curvatures via [@hornHillShadingReflectance1981], [@zevenbergen1987quantitative], and [@minarComprehensiveSystemDefinitions2020].
- **Relative terrain indices**: roughness, Topographic Position Index [@wilsonMultiscaleTerrainAnalysis2007], Bathymetric Position Index [@lundbladBenthicTerrainClassification2006], rugosity, entropy, percentile elevation [@lundquistAutomatedAlgorithmMapping2008], horizon angles, and Sky View Factor; at multiple scales.
- **Terrain classification**: Morphological Filters [@keqizhangProgressiveMorphologicalFilter2003, @pingelImprovedSimpleMorphological2013a], and skewness balancing [@bartelsDTMGenerationLIDAR2006; @bartelsThresholdfreeObjectGround2010].
- **Hydrological analysis**: depression filling [@barnesPriorityFloodOptimalDepressionFilling2014], flow accumulation (using D8 [@jensonExtractingTopographicStructure1988], D-infinity [@tarbotonNewMethodDetermination1997] or FD8 [@quinnPredictionHillslopeFlow1991] methods), and relative height models [@nobreHeightNearestDrainage2011].
- **Visualisation**: hillshading [@burroughPrinciplesGeographicalInformation2015], multi-directional hillshading [@mark1992multidirectional], and Perceptually Shaded Slope Maps [@pingelPerceptuallyShadedSlope2014a].
- **Cost-distance analysis**: friction-based spreading via Tomlin [@tomlin1983digital], Eastman [@eastman1989pushbroom], and fast sweeping [@zhaoFastSweepingMethod2005].

`Geomorphometry.jl` was designed to be used (and extended) by researchers studying geomorphometry on global scales. Furthermore, it was designed to also teach geomorphometry, by providing several algorithms for a task and an exhaustive, visual documentation.

# State of the field

Several mature tools exist for terrain analysis, but each involves trade-offs between performance and extensibility.
In Python, `xDEM` provides DEM analysis with a focus on elevation change detection and uncertainty, while `RichDEM` offers optimised hydrological algorithms.
In R, `MultiscaleDTM` [@ilichMultiscaleDTMOpensourcePackage2023] focuses on multiscale terrain metrics.
These Python and R packages are convenient but, being interpreted, can be slow on large datasets.
At the other end of the spectrum, PCRaster offers a comprehensive modelling language with hydrological functions in C++, while `WhiteboxTools` [@lindsayWhiteboxGATCase2016] provides a wide range of geospatial algorithms in Rust.
Both are fast, but the languages make them difficult to extend or modify for novel research by practitioners.
Indeed, this was the original trigger for the creation of this package.
Julia's design---expressive like Python, yet compiled to efficient native code---offers both qualities simultaneously, echoing the goals of the Julia manifesto [@bezansonJuliaFreshApproach2017].
<!-- `Geomorphometry.jl` consolidates functionality that is typically scattered---terrain filtering, terrain derivatives, hydrological routing, and cost-distance analysis---into a single, consistent, and extensible Julia package. -->

# Software design

`Geomorphometry.jl` operates on standard Julia `AbstractMatrix` types, making it compatible with any array-like data structure in the ecosystem.
Geospatial awareness is added through package extensions for `GeoArrays.jl` and `Rasters.jl` rather than hard dependencies, keeping the core lightweight.
The package is designed to be customizable and extendable, as shown in its core concepts:

- **Multiple algorithms.** Where several methods exist for the same task (e.g., slope, flow direction, cost-distance), users can select the most appropriate one. Algorithm selection uses Julia's multiple dispatch, making it straightforward to add new methods without modifying existing code---useful for both research and teaching.
- **Geospatial aware.** Package extensions automatically extract cell sizes and handle coordinate reference systems---including geographic (lat/lon) DEMs---preventing common user errors.
- **Multiscale analysis.** Relative terrain indices accept configurable stencil windows from `Stencils.jl`, enabling analysis at multiple spatial scales, on both CPU and GPU, extending the approach of @ilichMultiscaleDTMOpensourcePackage2023.

Terrain derivative implementations are validated against GDAL's `gdaldem`, and flow direction against `pyflwdir`, ensuring consistency with established tools.

# Research impact statement

The package was essential for analysis of DEMs on a global scale, and even creating new ones [@pronkDeltaDTMGlobalCoastal2024].
The hydrologic analysis methods were used to automatically detect levees in DEMs [@pronkAutomatedLeveeDetection2026].
Furthermore, the package was used to demonstrate geomorphometry in the Python, R and Julia languages in the software chapter in the upcoming second edition of the Geomorphometry book [@henglGeomorphometryConceptsSoftware2009].

# AI usage disclosure

AI tools are used in the creation of the software and documentation. Specifically, they are used to review implementations, find bugs or other inconsistencies. The same holds for documentation, where they are used for autocompletion, review and to fix grammar and spelling mistakes.
No AI tools were used for critical components, such as the overall design, architecture, or validation of the package.

# Acknowledgements

We thank the authors of the packages in Julia(Geo) ecosystem, particularly `GeoInterface`, `Stencils.jl`, and `DataStructures.jl`, on which Geomorphometry.jl is built.

# References
