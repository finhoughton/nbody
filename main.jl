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
        centre::Point,
        side_length::Float64
    )::Maybe{BHTree}
        
        len = length(particles)

        if isempty(particles)
            return EmptyMaybe
        end

        centre_of_mass::Point = sum(p -> p.pos, particles) / len
        half_side_len::Float64 = 0.5 * side_length
        quarter_side_len::Float64 = 0.5 * half_side_len 

        nws::Vector{Particle} = []
        nes::Vector{Particle} = []
        sws::Vector{Particle} = []
        ses::Vector{Particle} = []
        quads = (nws, nes, sws, ses)

        sizehint!.(quads, 1 + len ÷ 4)

        if length(particles) ≠ 1
            for particle in particles
                push_to_vector!(quads..., particle, centre)
            end
        end

        centres = (
            centre + Point(-quarter_side_len, quarter_side_len),
            centre + Point(quarter_side_len, quarter_side_len), 
            centre + Point(-quarter_side_len, -quarter_side_len),
            centre + Point(quarter_side_len, -quarter_side_len))

        subtrees = [Maybe(BHTree(quad, quad_centre, half_side_len)) for (quad, quad_centre) in zip(quads, centres)]
        new(particles, subtrees..., centre_of_mass, side_length)
    end
end

function push_to_vector!(
    nws::Vector{Particle},
    nes::Vector{Particle},
    sws::Vector{Particle},
    ses::Vector{Particle},
    p::Particle,
    c::Particle
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