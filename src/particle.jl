using StaticArrays
using LinearAlgebra

include("utils.jl")

mutable struct Particle
    mass::Float64
    pos::SVector{2, Float64}
    v::SVector{2, Float64}
    force_applied::SVector{2, Float64}
    fixed::Bool
end

function Particle(
    ;mass::Float64,
    pos::SVector{2, Float64}=zeros(Float64, 2),
    v::SVector{2, Float64}=zeros(Float64, 2),
    fixed::Bool=false
)::Particle
    if mass ≤ 0
        error("particle mass must be positive")
    elseif fixed && norm(v) ≠ 0
        error("fixed is incomaptable with velocity.")
    end
    Particle(mass, pos, v, zeros(Float64, 2), fixed)
end

function calculate_force(p::Particle, q::Particle)::SVector{2, Float64}
    normalize(q.pos - p.pos) * G * p.mass * (q.mass * inv(norm(p.pos - q.pos) ^ 2))
end

# ----- creaing particles -----

function random_particle() :: Particle
    mass::Float64 = M_EARTH * (rand() + 0.5)
    position::SVector{2, Float64} = SA[rand() - 0.5, rand() - 0.5] * DIST * 4
    velocity::SVector{2, Float64} = SA[rand() - 0.5, rand() - 0.5] * DIST * 0.2
    Particle(mass=mass, pos=position, v=velocity)
end

function random_particles(
    ;n::Int64,
    edge_len::Float64,
    mass_mean::Float64,
    mass_stddev::Float64,
    velocity_mean::Float64,
    velocity_stddev::Float64
    ) :: Vector{Particle}

    masses = Vector{Float64}(undef, n)
    positions = Vector{SVector{2, Float64}}(undef, n)
    velocities = Vector{SVector{2, Float64}}(undef, n)

    randn!(Float64, masses)
    map!(x -> abs(x) * mass_stddev + mass_mean, masses)

    rand!(SVector{2, Float64}, positions)
    map!(x -> x * edge_len, positions)

    randn!(SVector{2, Float64}, velocities)
    map!(x -> x * velocity_stddev .+ velocity_mean, velocities)

    particles::Vector{Particle} = Particle.(masses, positions, velocities)
end