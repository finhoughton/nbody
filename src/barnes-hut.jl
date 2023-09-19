using StaticArrays

include("particle.jl")

struct BHTree
    particles::Vector{Particle}
    NW::Maybe{BHTree}
    NE::Maybe{BHTree}
    SW::Maybe{BHTree}
    SE::Maybe{BHTree}
    centre_of_mass::SVector{2, Float64}
    side_length::Float64

    function BHTree(
        particles::Vector{Particle},
        centre::SVector{2, Float64},
        side_length::Float64
    )::Maybe{BHTree}
        
        len = length(particles)

        if len == 0
            return nothing
        end

        half_side_len::Float64 = 0.5 * side_length
        quarter_side_len::Float64 = 0.5 * half_side_len 

        quadrants::Tuple{Vector{Particle}, Vector{Particle}, Vector{Particle}, Vector{Particle}} = ([], [], [], [])

        if length(particles) ≠ 1
            total_mass::Float64 = 0
            total_pos::SVector{2, Float64} = SVector(0, 0)
            for particle ∈ particles
                total_mass += particle.mass
                total_pos += particle.pos
                push_to_vector!(quadrants..., particle, centre)
            end
            centre_of_mass::SVector{2, Float64} = total_pos / total_mass
        else
            centre_of_mass = first(particles).pos
        end

        centres = (
            centre + SVector(-quarter_side_len, quarter_side_len),
            centre + SVector(quarter_side_len, quarter_side_len), 
            centre + SVector(-quarter_side_len, -quarter_side_len),
            centre + SVector(quarter_side_len, -quarter_side_len))

        children = [BHTree(quad, quad_centre, half_side_len) for (quad, quad_centre) ∈ zip(quadrants, centres)]
        Just(new(particles, children..., centre_of_mass, side_length))
    end
end

function push_to_vector!(
    nws::Vector{Particle},
    nes::Vector{Particle},
    sws::Vector{Particle},
    ses::Vector{Particle},
    p::Particle,
    c::SVector{2, Float64}
)::Nothing
    if (p.pos[1] >= c[1] && p.pos[2] >= c[2])
        push!(nws, p)
    elseif (p.pos[1] < c[1] && p.pos[2] >= c[2])
        push!(nes, p)
    elseif (p.pos[1] >= c[1] && p.pos[2] < c[2])
        push!(sws, p)
    elseif (p.pos[1] < c[1] && p.pos[2] < c[2]) 
        push!(ses, p)
    end
    nothing
end