
include("maybe.jl")

struct Point
    x::Float64
    y::Float64
end

function +(a::Point, b::Point)::Point
    Point(a.x + b.x, a.y + b.y) 
end

function compare(a::Point, b::Point)::(Int, Int)
    cmp_x = a.x > b.x ? 1 : -1
    cmp_y = a.y > b.y ? 1 : -1
    (cmp_x, cmp_y)
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
        if isempty(particles)
            return Maybe{nothing}
        end

        centre_of_mass::Point = sum(map(p -> p.pos, particles)) / length(particles)
        half_side_len::Float64 = side_length / 2
        quarter_side_len::Float64 = half_side_len / 2

        if length(particles) ≠ 1
            nws::Vector{Particle} = filter(x -> compare(x, centre) == (1, 1), particles)
            nes::Vector{Particle} = filter(x -> compare(x, centre) == (-1, 1), particles)
            sws::Vector{Particle} = filter(x -> compare(x, centre) == (1, -1), particles)
            ses::Vector{Particle} = filter(x -> compare(x, centre) == (-1, -1), particles)
        else
            nws = nes = sws = ses = []
        end
        NW::BHTree = BHTree(nws, (centre + Point(-quarter_side_len, quarter_side_len)), half_side_len)
        NE::BHTree = BHTree(nes, (centre + Point(quarter_side_len, quarter_side_len)), half_side_len)
        SW::BHTree = BHTree(sws, (centre + Point(-quarter_side_len, -quarter_side_len)), half_side_len)
        SE::BHTree = BHTree(ses, (centre + Point(quarter_side_len, -quarter_side_len)), half_side_len)

        self::BHTree = new(particles, NW, NE, SW, SE, centre_of_mass, side_length)
        return Maybe(self)
    end
end
