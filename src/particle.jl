using StaticArrays
using LinearAlgebra

include("utils/maybe.jl")

mutable struct Particle
    mass::Float64
    pos::SVector{2, Float64}
    v::SVector{2, Float64}
    force_applied::SVector{2, Float64}
    fixed::Bool
    function Particle(
        mass::Float64,
        pos::SVector{2, Float64},
        v::SVector{2, Float64}=zeros(Float64, 2); 
        fixed::Bool=false
    )::Particle
        if size(v, 1) ≠ 2
            error("Particle only supports 2D positions and velocity")
        elseif mass ≤ 0
            error("particle mass must be positive")
        elseif fixed && norm(v) ≠ 0
            error("fixed is incomaptable with velocity.")
        end
        new(mass, pos, v, zeros(Float64, 2), fixed)
    end
end

