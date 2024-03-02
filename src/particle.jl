using StaticArrays
using LinearAlgebra

include("utils.jl")

mutable struct Particle
    id::Int
    mass::Float64
    pos::SVector{2,Float64}
    v::SVector{2,Float64}
    force_applied::SVector{2,Float64}
    fixed::Bool
end

function Particle(
    ; mass::Float64,
    pos::SVector{2,Float64}=SA[0.0, 0.0],
    v::SVector{2,Float64}=SA[0.0, 0.0],
    fixed::Bool=false
)::Particle
    if mass ≤ 0
        error("particle mass must be positive")
    elseif fixed && norm(v) ≠ 0
        error("fixed is incomaptable with velocity.")
    end
    Particle(-1, mass, pos, v, SA[0.0, 0.0], fixed)
end

Particle(mass::Float64, pos::SVector{2,Float64}, v::SVector{2,Float64}) = Particle(mass=mass, pos=pos, v=v, fixed=false)

function calculate_force(p::Particle, q::Particle)::SVector{2,Float64}
    normalize(q.pos - p.pos) * G * p.mass * (q.mass * inv(EPS_SOFTENING + norm(p.pos - q.pos)^2))
end

function Base.:(==)(p::Particle, q::Particle)::Bool
    p.id == q.id
end

# Compute net acceleration on particle p at position pos, with all other particles
# held at their current positions. Used for RK4 sub-step force evaluations.
function acceleration_at(p::Particle, particles::Vector{Particle}, pos::SVector{2,Float64})::SVector{2,Float64}
    if p.fixed
        return SA[0.0, 0.0]
    end
    a = SA[0.0, 0.0]
    for q ∈ particles
        if q == p
            continue
        end
        diff = q.pos - pos
        a += normalize(diff) * G * q.mass * inv(EPS_SOFTENING + norm(diff)^2)
        # mass of p cancels: F/m = G * q.mass / (eps + d^2)
    end
    return a
end

# Advance particle p by one timestep Δt using 4th-order Runge-Kutta integration.
# State vector is (pos, v); derivative is (v, a) where a = sum of gravitational
# accelerations from all other particles.
function rk4_update_particle!(max_speed::Float64, particles::Vector{Particle}, p::Particle)::Nothing
    if p.fixed
        p.force_applied = SA[0.0, 0.0]
        return nothing
    end

    pos0 = p.pos
    v0   = p.v

    # k1 — derivative at current state
    k1_v   = acceleration_at(p, particles, pos0)
    k1_pos = v0

    # k2 — derivative at midpoint estimated with k1
    k2_v   = acceleration_at(p, particles, pos0 + (Δt/2) * k1_pos)
    k2_pos = v0 + (Δt/2) * k1_v

    # k3 — derivative at midpoint estimated with k2
    k3_v   = acceleration_at(p, particles, pos0 + (Δt/2) * k2_pos)
    k3_pos = v0 + (Δt/2) * k2_v

    # k4 — derivative at end of interval estimated with k3
    k4_v   = acceleration_at(p, particles, pos0 + Δt * k3_pos)
    k4_pos = v0 + Δt * k3_v

    new_pos = pos0 + (Δt/6) * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)
    new_v   = v0   + (Δt/6) * (k1_v   + 2*k2_v   + 2*k3_v   + k4_v)

    if norm(new_v) > max_speed
        new_v = normalize(new_v) * max_speed
    end

    p.pos = new_pos
    p.v   = new_v
    p.force_applied = SA[0.0, 0.0]
    nothing
end

# this function is described in the pseudocode section of the design
function random_particles(
    ; n::Int64,
    edge_len::Float64,
    mass_mean::Float64,
    mass_stddev::Float64,
    velocity_stddev::Float64
)::Vector{Particle}

    masses::Vector{Float64} = map(x -> abs(x) * mass_stddev + mass_mean + 1, randn(n))
    # take absolute value because can't have negative mass, add 1kg to avoid 0 mass

    x_pos::Vector{Float64} = randn(n)
    y_pos::Vector{Float64} = randn(n)
    positions::Vector{SVector{2,Float64}} = map(x -> SVector{2,Float64}(x) * edge_len, zip(x_pos, y_pos))

    x_vel::Vector{Float64} = randn(n)
    y_vel::Vector{Float64} = randn(n)
    velocities::Vector{SVector{2,Float64}} = map(x -> SVector{2,Float64}(x) * velocity_stddev, zip(x_vel, y_vel))

    return Particle.(masses, positions, velocities)
end

