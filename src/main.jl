using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Plots
gr()

include("barnes-hut.jl")

# order of include():
# main
# barnes-hut
# particle
# utils

# ----- constants -----

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 3

const EPS_SOFTENING::Float64 = 1e10
# stop forcess becoming too big when objects are very close

const EDGE::Float64 = 10 * DIST

const X_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)
const Y_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)

# ------ displaying ------

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

function step!(particles::Vector{Particle}, root::BHTree, Δt::Float64)::Nothing
    the_q = Queue{BHTree}()
    for p ∈ particles
        step_particle!(p, root, the_q)
        a::SVector{2, Float64} = p.force_applied / p.mass
        dv::SVector{2, Float64} = a * Δt
        p.v += dv
        p.pos += p.v * Δt
    end
end

# ----- saving ------

const delimiter::String = "   "

function save_simulation!(file::IOStream, ps::Vector{Particle}, iteration_num::Integer = 0)::Nothing
    write(file, string(iteration_num, delimiter, length(ps), "\n"))
    fields::Tuple{Vararg{Symbol}} = fieldnames(Particle)
    for p ∈ ps
        join(file, [getfield(p, f) for f ∈ fields], delimiter)
        write(file, "\n")
    end
    nothing
end

function parse_particle(s::String)::Particle
    p = [type(eval(Meta.parse(data))) for (type, data) ∈ zip(Particle.types, split(s, delimiter))]
    Particle(p...)
end

function read_simulation(file::IOStream)::Tuple{Int64, Vector{Particle}}
    iteration_num, num_particles = parse.(Int64, tuple(string.(split(readline(file), delimiter))...))
    ps::Vector{Particle} = Vector{Particle}(undef, num_particles)
    for (idx, line) ∈ enumerate(eachline(file))
        ps[idx] =  parse_particle(line)
    end
    (iteration_num, ps)
end

# ----- main function -----

function random_particle() :: Particle
    mass::Float64 = M_EARTH * (rand() + 0.5)
    position = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 4
    velocity = SVector{2, Float64}(rand() - 0.5, rand() - 0.5) * DIST * 0.2
    Particle(mass=mass, pos=position, v=velocity)
end

function main()::Nothing
    particles::Vector{Particle} = [random_particle() for _ ∈ 1:5]
    t::Float64 = 0
    Δt::Float64 = 1/32
    f = open("save.txt", "w")
    save_simulation!(f, particles)
    close(f)
    f = open("save.txt", "r")
    ps = read_simulation(f)
    close(f)
    # while !isempty(particles)
    #     start::DateTime = now()
    #     root = unsafe_from_just(BHTree(particles, SVector{2, Float64}(0, 0), 2 * EDGE))
    #     step!(particles, root, Δt)
    #     filter!(p -> maximum(abs, p.pos) < EDGE, particles) # delete offscreen particles
    #     showparticles(particles)

    #     timetaken = convert(Millisecond, now() - start)
    #     if (v = value(timetaken)) < Δt * 1000
    #         sleeptime = Δt - v/1000
    #         println("sleeping for $sleeptime")
    #         sleep(sleeptime)
    #     end
    #     t += Δt
    # end
    nothing
end

main()
