module GeomorphometryEikonalExt

using Geomorphometry
using Eikonal

"""
    spread(::FastSweeping, points::AbstractMatrix{<:Real}, initial::AbstractMatrix{<:Real}, friction::AbstractMatrix{<:Real}; )

Total friction distance spread from `points` as described by [Zhao, H (2005)](@cite zhaoFastSweepingMethod2005). 

# Output
- `Array{Float64,2}` Total friction distance

# Arguments
- `points::Vector{CartesianIndex}` Input Array
- `initial::AbstractVector{<:Real}` Initial values of the result
- `friction::Matrix{<:Real}` Friction map
"""
function Geomorphometry.spread(
    fs::Geomorphometry.FastSweeping,
    points,
    initial,
    friction;
    kwargs...,
)
    result = similar(friction)

    solver = Eikonal.FastSweeping(parent(friction))
    for I in points
        solver.t[I] = initial[I]
    end
    Eikonal.sweep!(solver; nsweeps = fs.iterations, verbose = fs.debug, epsilon = fs.eps)
    result .= solver.t[2:end, 2:end]
    result
end

end
