# Cloth Simulation Filter (Zhang et al., 2016)"
# using Plots

const gravity = 9.85
const smoothThreshold = 0.3
CI = CartesianIndex

const nb = (
    (CI(0, 0), CI(1, 0)),
    (CI(0, 0), CI(0, 1)),
    (CI(0, 0), CI(1, 1)),
    (CI(0, 0), CI(-1, 1)),
    (CI(0, 0), CI(1, -1)),
    (CI(0, 0), CI(2, 0)),
    (CI(0, 0), CI(0, 2)),
    (CI(0, 0), CI(2, 2)),
    (CI(0, 0), CI(-2, 2)),
    (CI(0, 0), CI(2, -2)),
)

function csf(A; iterations = 500, time_step = 0.25, rigidness = 3, postprocess = true, threshold = 0.05, damping = 0.01)

    min_A, max_A = extrema(skipmissing(A))
    cloth = copy(A)

    nodata = ismissing.(A)
    moveable = .~nodata
    cloth[moveable] .= min_A
    prev_cloth = copy(cloth)
    acceleration = gravity * time_step * time_step
    temp = zero(eltype(A))

    rigidness_factor = (1 - 0.5^rigidness)

    R = CartesianIndices(cloth)
    w, h = size(A)

    max_diff = zero(eltype(A))
    prev_max_diff = typemax(max_diff)

    @showprogress 1 "Cloth simulating..." for i ∈ 1:iterations

        # Drop cloth each iteration (gravity)
        @inbounds for I in R
            if moveable[I]
                temp = cloth[I]
                move = (cloth[I] - prev_cloth[I]) * (1.0 - damping) + acceleration * time_step * time_step
                cloth[I] += move
                prev_cloth[I] = temp
                # println("$move $(cloth[I] - prev_cloth[I]) $(acceleration * time_step * time_step)")
            end
        end
        # plot((
        #         heatmap(A; colorbar = false, clims = (min_A, max_A), c = :viridis),
        #         heatmap(cloth; colorbar = false, clims = (min_A, max_A), c = :viridis),
        #         heatmap(moveable; colorbar = false, clims = (0, 1), c = :grays)
        #     )...,
        #     layout = (1, 3), ratio = :equal)


        # If we've reached the ground, set node to unmoveable
        max_diff = zero(eltype(A))
        @inbounds for I in R
            if moveable[I]
                diff = abs(prev_cloth[I] - cloth[I])
                max_diff = max(max_diff, diff)

                if cloth[I] > A[I]
                    cloth[I] = A[I]
                    moveable[I] = false
                end
            end
        end
        # println("$max_diff $prev_max_diff")

        # Stop if the cloth stopped moving
        if (max_diff < (threshold / 100)) || (abs(prev_max_diff - max_diff) < (threshold / 100))
            @info "Done in $i iterations."
            break
        end
        prev_max_diff = max_diff

        # Move each cloth node based on spring interaction with neighbors
        @inbounds for I in CartesianIndices((3:w-2, 3:h-2))

            # Only get graph, so
            for (a, b) in nb
                ai, bi = I + a, I + b

                # vertical difference between nodes
                correction = cloth[bi] - cloth[ai]
                correction *= rigidness_factor

                if (moveable[ai] && moveable[bi])
                    cloth[ai] += correction * 0.5
                    cloth[bi] -= correction * 0.5
                elseif (moveable[ai] && !moveable[bi] && !nodata[bi])
                    cloth[ai] += correction
                elseif (!moveable[ai] && !nodata[ai] && moveable[bi])
                    cloth[bi] -= correction
                end
            end
        end



        # plot((
        #         heatmap(A; colorbar = false, clims = (min_A, max_A), c = :viridis),
        #         heatmap(cloth; colorbar = false, clims = (min_A, max_A), c = :viridis),
        #         heatmap(moveable; colorbar = false, clims = (0, 1), c = :grays)
        #     )...,
        #     layout = (1, 3), ratio = :equal)

        # if (params.bSloopSmooth) {
        #     cloth.movableFilter();
        # }
    end

    return cloth, moveable
end
