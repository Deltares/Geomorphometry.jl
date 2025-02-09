module GeomorphometryGeoArraysExt

using Geomorphometry
using GeoArrays
using FillArrays

function Geomorphometry.xyratio(dem::GeoArray)
    return Fill(1.0, size(dem))
end

end # module
