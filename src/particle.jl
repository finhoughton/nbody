using StaticArrays
using LinearAlgebra

include("utils.jl")

mutable struct Particle
    id::Int
    mass::Float64
    pos::SVector{2, Float64}
    v::SVector{2, Float64}
    previous_force_applied::SVector{2, Float64}
    force_applied::SVector{2, Float64}
    fixed::Bool
end


function Particle(
    ;mass::Float64,
    pos::SVector{2, Float64}=SA[0.0, 0.0],
    v::SVector{2, Float64}=SA[0.0, 0.0],
    fixed::Bool=false
)::Particle
    if mass ≤ 0
        error("particle mass must be positive")
    elseif fixed && norm(v) ≠ 0
        error("fixed is incomaptable with velocity.")
    end
    Particle(-1, mass, pos, v, SA[0.0, 0.0], SA[0.0, 0.0], fixed)
end

Particle(mass::Float64, pos::SVector{2, Float64}, v::SVector{2, Float64}) = Particle(mass=mass, pos=pos, v=v, fixed=false)

function calculate_force(p::Particle, q::Particle)::SVector{2, Float64}
    normalize(q.pos - p.pos) * G * p.mass * (q.mass * inv(EPS_SOFTENING + norm(p.pos - q.pos) ^ 2))
end

function Base.:(==)(p::Particle, q::Particle)::Bool
    p.id == q.id
end

# iterate the particle's velocity and position
function update_particle!(p::Particle)::Nothing
    # area of a trapezium = (a + b) * h / 2
    if !p.fixed
        # using trapezium approximation for impulse (I). I = ∫F dt
        impulse::SVector{2, Float64} = 0.5 * Δt * (p.previous_force_applied + p.force_applied)
        Δv::SVector{2, Float64} = impulse / p.mass 

        old_v = p.v
        p.v += Δv
        new_v = p.v

        # use trapezium approximation for displacement (s). s = ∫v dt
        Δs = 0.5 * Δt * (old_v + new_v)
        p.pos += Δs
    end 
    p.previous_force_applied = p.force_applied
    p.force_applied = SA[0.0, 0.0]
    nothing
end

# this function is described in the pseudocode section of the design
function random_particles(
    ;n::Int64,
    edge_len::Float64,
    mass_mean::Float64,
    mass_stddev::Float64,
    velocity_mean::Float64,
    velocity_stddev::Float64
    ) :: Vector{Particle}

    masses::Vector{Float64} = map(x -> abs(x) * mass_stddev + mass_mean + 1, randn(n))
    # take absolute value because can't have negative mass, add 1kg to avoid 0 mass

    x_pos::Vector{Float64} = randn(n)
    y_pos::Vector{Float64} = randn(n)
    positions::Vector{SVector{2, Float64}} = map(x -> SVector{2, Float64}(x) * edge_len, zip(x_pos, y_pos))

    x_vel::Vector{Float64} = randn(n)
    y_vel::Vector{Float64} = randn(n)
    velocities::Vector{SVector{2, Float64}} = map(x -> SVector{2, Float64}(x) * velocity_stddev .+ (velocity_mean / sqrt(2)), zip(x_vel, y_vel))

    return Particle.(masses, positions, velocities)
end

