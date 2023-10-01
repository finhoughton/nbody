using LinearAlgebra
using StaticArrays

include("barnes-hut.jl")

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 30

const EPS_SOFTENING::Float64 = 1e7
# stop forcess becoming too big when objects are very close

const EDGE::Float64 = 10 * DIST

const X_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)
const Y_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)

function test1(p::Vector{Particle})
    root = BHTree(p, SVector{2, Float64}(0, 0), 1e10)
    nothing
end

function test2(p::Vector{Particle})
    root = BHTree(1, p, SVector{2, Float64}(0, 0), 1e10)
    nothing
end

function random_particle() :: Particle
    m::Float64 = M_EARTH * (rand() + 0.5)
    pos = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 4
    v = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 0.5
    Particle(m, pos, v)
end
 
function main()::Nothing
    particles::Vector{Particle} = [random_particle() for _ ∈ 1:100000]
    root = unsafe_from_just(BHTree(particles, SVector{2, Float64}(0, 0), 1e10))
    step!(particles, root)
    print("ran")
    nothing
end

main()
