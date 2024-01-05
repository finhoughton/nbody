using DataStructures

include("particle.jl")

const MAC::Float64 = 1 # smaller MAC = more accurate

struct BHTree
    particles::Vector{Particle}
    children::SVector{4, Maybe{BHTree}}
    NW::Maybe{BHTree}
    NE::Maybe{BHTree}
    SW::Maybe{BHTree}
    SE::Maybe{BHTree}
    centre::SVector{2, Float64}
    total_mass::Float64
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

        quadrants::Tuple{Vector{Particle}, Vector{Particle}, Vector{Particle}, Vector{Particle}} = ([], [], [], [])

        if len == 1
            total_mass = first(particles).mass
            centre_of_mass = first(particles).pos
        else
            # centre of mass = (m_1r_1 + m_2r_2 + ...)/(m_1 + m_2 + ...)
            total_mass::Float64 = 0
            total::SVector{2, Float64} = SA{Float64}[0.0, 0.0]
            for particle ∈ particles
                total_mass += particle.mass
                total += particle.pos * particle.mass
                push_to_quadrant!(quadrants..., particle, centre)
            end
            centre_of_mass::SVector{2, Float64} = total / total_mass
        end

        half_side_len::Float64 = 0.5 * side_length
        quarter_side_len::Float64 = 0.5 * half_side_len 
        centres = (
            centre + SVector(-quarter_side_len, quarter_side_len),  # NW
            centre + SVector(quarter_side_len, quarter_side_len),   # NE
            centre + SVector(-quarter_side_len, -quarter_side_len), # SW
            centre + SVector(quarter_side_len, -quarter_side_len))  # SE

        # recurive call creating the 4 children
        children = SVector{4, Maybe{BHTree}}([BHTree(quad, quad_centre, half_side_len) for (quad, quad_centre) ∈ zip(quadrants, centres)])
        return Just(new(particles, children, children..., centre, total_mass, centre_of_mass, side_length))
    end
end

function push_to_quadrant!(
    nws::Vector{Particle},
    nes::Vector{Particle},
    sws::Vector{Particle},
    ses::Vector{Particle},
    p::Particle,
    c::SVector{2, Float64}
    )::Nothing

    if p.pos[1] >= c[1]
        # it is west of the centre
        if p.pos[2] >= c[2]
            push!(nws, p)
        else
            push!(sws, p)
        end

    else
        # it is east of the centre
        if p.pos[2] >= c[2]
            push!(nes, p)
        else
            push!(ses, p)
        end
    end
    nothing
end


function calculate_force(p::Particle, node::BHTree)::SVector{2, Float64}
    normalize(node.centre_of_mass - p.pos) * G * p.mass * (node.total_mass * inv(EPS_SOFTENING + norm(p.pos - node.centre_of_mass) ^ 2))
end
