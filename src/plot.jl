# using RecipesBase
using ImageFiltering
using ImageCore

"""
```
image = pssm(A; exaggeration, resolution)
```
Perceptually Shaded Slope Map by *Pingel, Clarke. 2014* [^pingel2014].

# Output
- `image::Gray{T,2}` Grayscale image

# Arguments
- `A::Array{Real,2}` Input Array
- `exaggeration::Real=2.3` Factor to exaggerate elevation
- `resolution::Real=1.0` Resolution of cell size

[^pingel2014]: Pingel, Thomas, and Clarke, Keith. 2014. ‘Perceptually Shaded Slope Maps for the Visualization of Digital Surface Models’. Cartographica: The International Journal for Geographic Information and Geovisualization 49 (4): 225–40. <https://doi.org/10/ggnthv>.
"""
function pssm(A::AbstractMatrix{<:Real}; exaggeration=2.3, resolution=1.)
    x, y = imgradients(A * exaggeration, ImageFiltering.sobel)
    G = sqrt.(x.^2 .+ y.^2)

    G /= resolution  # account for horizontal resolution

    f = scaleminmax(0, 1)
    clamped = f.(G)
    Gray.(1 .- clamped)
end
