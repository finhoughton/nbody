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