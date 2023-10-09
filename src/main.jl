using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Plots
gr()

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

function random_particle() :: Particle
    m::Float64 = M_EARTH * (rand() + 0.5)
    pos = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 4
    v = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 0.5
    Particle(m, pos, v)
end

function showparticles(particles::Vector{Particle})::Nothing

    positions::Vector{Vector{Float64}} = [p.pos for p in particles]

    xs::Vector{Float64} = getindex.(positions, 1)
    ys::Vector{Float64} = getindex.(positions, 2)
    p::Plots.Plot = scatter(xs, ys, legend=false, showaxis=false, background_colour=:black)
    xlims!(p, X_LIMITS...)
    ylims!(p, Y_LIMITS...)
    display(p)
    nothing
end

function step_particle!(p::Particle, root::BHTree, the_q::Queue{BHTree}):: Nothing
    enqueue!(the_q, root)
    while !isempty(the_q)
        current::BHTree = dequeue!(the_q)
        distance_to_centre::Float64 = norm(p.pos - current.centre)
        if current.side_length < MAC * distance_to_centre
            p.force_applied += calculate_force(p, current)
        else
            from_maybe_with.(
                nothing,
                x -> enqueue!(the_q, x),
                current.children
                )
        end
    end
    nothing
end

function step!(particles::Vector{Particle}, root::BHTree, dt::Float64)::Nothing
    the_q = Queue{BHTree}()
    for p ∈ particles
        step_particle!(p, root, the_q)
        a::SVector{2, Float64} = p.force_applied / p.mass
        dv::SVector{2, Float64} = a * dt
        p.v += dv
        p.pos += p.v * dt
    end
end

function main()::Nothing
    particles::Vector{Particle} = [random_particle() for _ ∈ 1:100]
    root = unsafe_from_just(BHTree(particles, SVector{2, Float64}(0, 0), EDGE))
    t::Float64 = 0
    dt::Float64 = 1/32
    while true
        start::DateTime = now()
        root = unsafe_from_just(BHTree(particles, SVector{2, Float64}(0, 0), EDGE))
        step!(particles, root)
        filter!(p -> maximum(abs, p.pos) < EDGE, particles) # delete offscreen particles
        showparticles(particles)

        timetaken = convert(Millisecond, now() - start)
        if (v = value(timetaken)) < dt * 1000
            a = dt - v/1000
            println("sleeping for $a")
            sleep(a)
        end
        t += dt
    end
    nothing
end

main()
