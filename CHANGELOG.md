# Changelog

All notable changes to this project will be documented in this file.

## [0.7.0] - 2025-03-21

### Added
- Improved documentation by using VitePress and showcasing all methods on a DEM.
- Added hydrology operators: `priorityflood`, `streamflow`
- Added *multiscale* options to some filters using a `window` kwarg for a Stencil from Stencils.jl package. Other methods now take a `radius` kwarg.
- Added `BPI`
- Overhauled curvature methods by introducing `plan_curvature`, `profile_curvature` and `contour_curvature`, while deprecating `curvature` for `laplacian`
- Added direction kwarg to slope and curvature methods.
- Added this `CHANGELOG.md`
- Added package extensions on GeoArrays, Rasters to support automatic cellsizes.
- Added package extension on Eikonal for a faster `spread` method.

### Changed
- Added Stencils as dependency.
- Refactored spread to choose from multiple algorithms

### Fixed
- Relaxed `Array` input to `AbstractArray` for `opening`
## [0.6.0] - 2023-08-25
### Changed
- Renamed package to Geomorphometry

### Added
- Added `erosion` parameter in `PMF`

## [0.5.2] - 2023-08-14
### Changed
- Light maintenance to CI scripts
- Compat updates

## [0.5.1] - 2023-01-12
### Added
- Relaxed input for `PMF`

## [0.5.0] - 2022-11-02

### Added
- Added `multihillshade`

### Changed
- Uses LocalFilters for terrain filters, improving performance, but changing the edge behaviour of some filters.

## [0.4.0] - 2022-10-17
### Added
- Added `hillshade`, `curvature`

### Changed
- Used LocalFilters for `PMF`.

## [0.3.2] - 2022-09-07
### Added
- Added terrain kernels for different algorithms.

### Fixed
- Fixed bug in TPI

## [0.3.1] - 2022-06-15
### 
- Added `skewness` filter

## [0.3.0] - 2022-02-01
### Added
- Added PSF (Progressive Slope Filter) function.

### Changed
- Improved performance of mapwindow function used in most filters.

## [0.2.0] - 2021-11-22
### Added
- Initial release of the GeoArrayOps package.
- Basic terrain analysis functions.

