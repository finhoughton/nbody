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
    pos::SVector{2, Float64}=SVector{2, Float64}(zeros(Float64, 2)),
    v::SVector{2, Float64}=SVector{2, Float64}(zeros(Float64, 2)),
    fixed::Bool=false
)::Particle
    if mass ≤ 0
        error("particle mass must be positive")
    elseif fixed norm(v) ≠ 0
        error("fixed is incomaptable with velocity.")
    end
    Particle(mass, pos, v, SVector{2, Float64}(zeros(Float64, 2)), fixed)
end

function calculate_force(p::Particle, q::Particle)::SVector{2, Float64}
    normalize(q.pos - p.pos) * G * p.mass * (q.mass * inv(EPS_SOFTENING + norm(p.pos - q.pos) ^ 2))
end

# ----- creaing particles -----

function random_particle() :: Particle
    mass::Float64 = M_EARTH * (rand() + 0.5)
    position::SVector{2, Float64} = SA[rand() - 0.5, rand() - 0.5] * DIST * 4
    velocity::SVector{2, Float64} = SA[rand() - 0.5, rand() - 0.5] * DIST * 2
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

    masses = randn(n)
    map!(x -> abs(x) * mass_stddev + mass_mean, masses)

    x_pos = randn(n)
    y_pos = randn(n)
    positions::Vector{SVector{2, Float64}} = map(SVector{2, Float64}, zip(x_pos, y_pos))
    map!(x -> x * edge_len, positions)

    x_vel = randn(n)
    y_vel = randn(n)
    velocities::Vector{SVector{2, Float64}} = map(SVector{2, Float64}, zip(x_vel, y_vel))
    map!(x -> x * velocity_stddev .+ velocity_mean, velocities)

    return Particle.(masses, positions, velocities)
end