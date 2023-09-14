using Match
using LinearAlgebra

include("utils/maybe.jl")

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 30

const EPS_SOFTENING::Float64 = 1e7
# stop forcess becoming too big when objects are very close

const EDGE::Float64 = 10 * DIST

const X_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)
const Y_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)

struct Point
    x::Float64
    y::Float64
end

function Base.:+(a::Point, b::Point)::Point
    Point(a.x + b.x, a.y + b.y) 
end

function Base.:-(a::Point, b::Point)::Point
    Point(a.x - b.x, a.y - b.y)
end

function Base.:*(a::Point, b::Real)::Point
    Point(a.x * b, a.y * b)
end

mutable struct Particle
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
        new(mass, pos, v, zeros(Float64, 2), fixed)
    end
end

function calculate_centre_of_mass(ps::Vector{Particle})::Point
    total_mass::Float64 = 0
    total_x::Float64 = 0
    total_y::Float64 = 0
    for p ∈ ps
        total_mass += p.mass
        total_x += p.pos.x
        total_y += p.pos.y
    end
    Point(total_x / total_mass, total_y / total_mass)
end

struct BHTree
    particles::Vector{Particle}
    NW::Maybe{BHTree}
    NE::Maybe{BHTree}
    SW::Maybe{BHTree}
    SE::Maybe{BHTree}
    centre_of_mass::Point
    side_length::Float64

    function BHTree(
        particles::Vector{Particle},
        centre::Point,
        side_length::Float64
    )::Maybe{BHTree}
        
        len = length(particles)

        if len == 0
            return nothing
        end

        centre_of_mass::Point = calculate_centre_of_mass(particles)
        half_side_len::Float64 = 0.5 * side_length
        quarter_side_len::Float64 = 0.5 * half_side_len 

        nws::Vector{Particle} = []
        nes::Vector{Particle} = []
        sws::Vector{Particle} = []
        ses::Vector{Particle} = []
        quads = (nws, nes, sws, ses)

        sizehint!.(quads, 3 + (len >> 2))

        if length(particles) ≠ 1
            for particle ∈ particles
                push_to_vector!(quads..., particle, centre)
            end
        end

        centres = (
            centre + Point(-quarter_side_len, quarter_side_len),
            centre + Point(quarter_side_len, quarter_side_len), 
            centre + Point(-quarter_side_len, -quarter_side_len),
            centre + Point(quarter_side_len, -quarter_side_len))

        children = [BHTree(quad, quad_centre, half_side_len) for (quad, quad_centre) in zip(quads, centres)]
        Just(new(particles, children..., centre_of_mass, side_length))
    end
end

function push_to_vector!(
    nws::Vector{Particle},
    nes::Vector{Particle},
    sws::Vector{Particle},
    ses::Vector{Particle},
    p::Particle,
    c::Point
)::Nothing
    if (p.pos.x >= c.x && p.pos.y >= c.y)
        push!(nws, p)
    elseif (p.pos.x < c.x && p.pos.y >= c.y)
        push!(nes, p)
    elseif (p.pos.x >= c.x && p.pos.y < c.y)
        push!(sws, p)
    elseif (p.pos.x < c.x && p.pos.y < c.y) 
        push!(ses, p)
    end
    nothing
end

function random_particle() :: Particle
    m::Float64 = M_EARTH * (rand() + 0.5)
    pos::Point = Point(rand() - 0.5, rand() - 0.5) * DIST * 4
    v::Vector{Float64} = (rand(Float64, 2) * 2 - ones(2)) * DIST * 0.5
    Particle(m, pos, v)
end

function main()::Nothing
    particles::Vector{Particle} = [random_particle() for _ ∈ 1:1000]
    root = BHTree(particles, Point(0, 0), 1e10)
    nothing
end

main()
