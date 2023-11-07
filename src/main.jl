using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Plots
using Random
gr()

include("barnes-hut.jl")

# order of include():
# main
# barnes-hut
# particle
# utils

# ----- constants -----

const FPS::Int64 = 50

const Δt::Float64 = 1/FPS

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 30

const EPS_SOFTENING::Float64 = 1e8
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

function step!(particles::Vector{Particle}, root::BHTree)::Nothing
    the_q = Queue{BHTree}()
    for p ∈ particles
        step_particle!(p, root, the_q)
        a::SVector{2, Float64} = p.force_applied / p.mass
        dv::SVector{2, Float64} = a * Δt
        p.v += dv
        p.pos += p.v * Δt
        p.force_applied = SA[0.0, 0.0]
    end
end

# ----- saving ------

const delimiter::String = "   "

function save!(file::IOStream, ps::Vector{T}, data::Vector{Int64})::Nothing where {T}
    push!(data, length(ps))
    join(file, data, delimiter)
    write(file, "\n")
    fields::Tuple{Vararg{Symbol}} = fieldnames(Particle)
    for p ∈ ps
        join(file, [getfield(p, f) for f ∈ fields], delimiter)
        write(file, "\n")
    end
    nothing
end

function read(file::IOStream, T::DataType)::Tuple{Tuple{Vararg{Int64}}, Vector{T}}
    (data..., num_particles) = parse.(Int64, tuple(string.(split(readline(file), delimiter))...))
    ps::Vector{T} = Vector{T}(undef, num_particles)
    for (idx, line) ∈ enumerate(eachline(file))
        particle = [type(eval(Meta.parse(data))) for (type, data) ∈ zip(T.types, split(line, delimiter))]
        ps[idx] = Particle(particle...)
    end
    (data, ps)
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

# ----- main function -----


function test_saving(particles::Vector{Particle})::Nothing
    fp = "data/save.txt"
    file = open(fp, "w")
    save!(file, particles, Vector{Int64}() )
    close(file)
    particles
    file = open(fp, "r")
    data, ps = read(file, Particle)
    ps::Vector{Particle}
    close(file)
    for (p1, p2) ∈ zip(ps, particles)
        println(p1)
        println(p2)
        println()
    end
    nothing
end

function main()::Nothing
    particles::Vector{Particle} = [random_particle() for _ ∈ 1:100]
    t::Float64 = 0
    while !isempty(particles)
        start::DateTime = now()
        root = unsafe_from_just(BHTree(particles, SA[0.0, 0.0], 2 * EDGE))
        step!(particles, root)
        filter!(p -> maximum(abs, p.pos) < EDGE, particles) # delete offscreen particles
        showparticles(particles)

        timetaken = convert(Millisecond, now() - start)
        if (v = value(timetaken)) < Δt * 1000
            sleeptime = Δt - v/1000
            println("sleeping for $sleeptime")
            sleep(sleeptime)
        end
        t += Δt
    end
    nothing
end

main()
