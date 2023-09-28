using DataStructures

include("particle.jl")

const MAC::Float64 = 1
const NOOP = () -> begin end

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

        half_side_len::Float64 = 0.5 * side_length
        quarter_side_len::Float64 = 0.5 * half_side_len 

        quadrants::Tuple{Vector{Particle}, Vector{Particle}, Vector{Particle}, Vector{Particle}} = ([], [], [], [])

        if length(particles) == 1
            total_mass = first(particles).mass
            centre_of_mass = first(particles).pos
        else
            # centre of mass = (m_1 * r_1 + m_2 * r_2 * ...)/(m_1 + m_2 + ...)
            total_mass::Float64 = 0
            total::SVector{2, Float64} = SVector(0, 0)
            for particle ∈ particles
                total_mass += particle.mass
                total += particle.pos * particle.mass
                push_to_vector!(quadrants..., particle, centre)
            end
            centre_of_mass::SVector{2, Float64} = total / total_mass
        end

        centres = (
            centre + SVector(-quarter_side_len, quarter_side_len),
            centre + SVector(quarter_side_len, quarter_side_len), 
            centre + SVector(-quarter_side_len, -quarter_side_len),
            centre + SVector(quarter_side_len, -quarter_side_len))

        children = SVector{4, Maybe{BHTree}}([BHTree(quad, quad_centre, half_side_len) for (quad, quad_centre) ∈ zip(quadrants, centres)])
        Just(new(particles, children, children..., centre, total_mass, centre_of_mass, side_length))
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

function calculate_force(p::Particle, node::BHTree)::SVector{2, Float64}
    normalize(node.centre_of_mass - p.pos) * G * p.mass * (node.total_mass * (1/norm(p.pos - node.centre_of_mass) ^ 2))
end

function calculate_force(p::Particle, q::Particle)::SVector{2, Float64}
    normalize(q.pos - p.pos) * G * p.mass * (q.mass * (1/norm(p.pos - q.pos) ^ 2))
end

function step_particle!(p::Particle, root::BHTree):: Nothing
    the_q = Queue{BHTree}()
    enqueue!(the_q, root)
    while length(the_q) ≠ 0
        current::BHTree = dequeue!(the_q)
        distance_to_centre::Float64 = norm(p.pos - current.centre)
        theta::Float64 = current.side_length / distance_to_centre
        if theta < MAC
            p.force_applied += calculate_force(p, current)
            continue
        end
        for c ∈ current.children
            if_just_then(NOOP, x -> enqueue!(the_q, x), c)
        end
    end
    nothing
end

function step!(particles::Vector{Particle}, root::BHTree)::Nothing
    for p ∈ particles
        step_particle!(p, root)
    end
    nothing
end
