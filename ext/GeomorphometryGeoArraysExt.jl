module GeomorphometryGeoArraysExt

using Geomorphometry
using GeoArrays
using FillArrays

degwidth::Float64 = 111_000.0

function Geomorphometry.cellsize(dem::GeoArray)
    T = GeoArrays.GI.crstrait(dem)
    _cellsize(T, dem)
end

_cellsize(::GeoArrays.GI.AbstractProjectedTrait, dem::GeoArray) = (dem.f.linear[1], dem.f.linear[4])
function _cellsize(::GeoArrays.GI.AbstractGeographicTrait, dem::GeoArray)
    centercoords = GeoArrays.coords(dem, round.(Int, size(dem)[1:2]./2))
    (dem.f.linear[1] * degwidth * cosd(centercoords[2]),
    dem.f.linear[4] * degwidth)
end

end # module
