using Match

include("maybe.jl")

struct Point
    x::Float64
    y::Float64
end

function +(a::Point, b::Point)::Point
    Point(a.x + b.x, a.y + b.y) 
end

let id::Int = 0
    mutable struct Particle
        id::Int
        mass::Float64
        pos::Point
        v::Vector{Float64}
        force_applied::Vector{Float64}
        fixed::Bool
        function Particle(
            mass::Float64,
            pos::Point,
            v::Vector{Float64}=zeros(Float64, 2);
            fixed::Bool=false
        )::Particle
            if size(v, 1) != 2
                error("Particle only supports 2D positions and velocity")
            elseif mass <= 0
                error("particle mass must be positive")
            elseif fixed && norm(v) != 0
                error("fixed is incomaptable with velocity.")
            end
            id += 1
            new(id, mass, pos, v, zeros(Float64, 2), fixed)
        end
    end
end

struct BHTree
    particles::Vector{Particle}
    NW::Maybe{BHTree}
    NE::Maybe{BHTree}
    SW::Maybe{BHTree}
    SE::Maybe{BHTree}
    centre_of_mass::Maybe{Vector{Float64}}
    side_length::Float64
    function BHTree(
        particles::Vector{Particle},
        c::Point,
        side_length::Float64
    )::Maybe{BHTree}

        if isempty(particles)
            return Maybe{nothing}
        end

        function push_to_vector!(
            nws::Vector{Particle},
            nes::Vector{Particle},
            sws::Vector{Particle},
            ses::Vector{Particle},
            p::Particle
        )::Nothing
            if (p.x >= c.x && p.y >= c.y)
                push!(nws, particle)
            elseif (p.x < c.x && p.y >= c.y)
                push!(nes, particle)
            elseif (p.x >= c.x && p.y < c.y)
                push!(sws, particle)
            elseif (p.x < c.x && p.y < c.y) 
                push!(ses, particle)
            end
        end

        centre_of_mass::Point = sum(p -> p.pos, particles) / length(particles)
        half_side_len::Float64 = side_length / 2
        quarter_side_len::Float64 = half_side_len / 2

        approx_size = length(particles) ÷ 4
        nws::Vector{Particle} = []
        nes::Vector{Particle} = []
        sws::Vector{Particle} = []
        ses::Vector{Particle} = []
        sizehint!(nws, approx_size)
        sizehint!(nes, approx_size)
        sizehint!(sws, approx_size)
        sizehint!(ses, approx_size)

        if length(particles) ≠ 1
            for particle in particles
                push_to_vector!(nws, nes, sws, ses, particle)
            end
        end

        NW::BHTree = Maybe(BHTree(nws, (centre + Point(-quarter_side_len, quarter_side_len)), half_side_len))
        NE::BHTree = Maybe(BHTree(nes, (centre + Point(quarter_side_len, quarter_side_len)), half_side_len))
        SW::BHTree = Maybe(BHTree(sws, (centre + Point(-quarter_side_len, -quarter_side_len)), half_side_len))
        SE::BHTree = Maybe(BHTree(ses, (centre + Point(quarter_side_len, -quarter_side_len)), half_side_len))

        new(particles, NW, NE, SW, SE, centre_of_mass, side_length)
    end
end
