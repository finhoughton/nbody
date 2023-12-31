using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Plots
using Combinatorics: combinations
using Random
gr()

include("iteration-algorithms.jl")

# order of include():
# main
# iteration algorithms
# particle
# utils

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

# ------ iteration ------

# barnes-hut
function step_particle!(root::BHTree, the_q::Queue{BHTree}, p::Particle)::Nothing 
    enqueue!(the_q, root)
    # enqueue the root node
    while !isempty(the_q)
        current::BHTree = dequeue!(the_q)
        # dequeue a node
        distance_to_centre::Float64 = norm(p.pos - current.centre)
        # calculate the distance to the centre to be used in MAC calculations
        if current.side_length < MAC * distance_to_centre
            # the ratio is not greater than the MAC, use this quadrant to appriximate the force on the particle
            p.force_applied += calculate_force(p, current)
        else
            # the ratio is greater than the MAC, enqueue the current node's children
            from_maybe_with.(
                nothing,
                x -> enqueue!(the_q, x),
                current.children
                )
            # enqueueing the children is a bit weird because they are `Maybe`s
        end
    end
    nothing
end

# barnes-hut step
function step!(particles::Vector{Particle}, root::BHTree)::Nothing
    the_q = Queue{BHTree}()
    # create the queue used in bh algorithm
    foreach(partial(step_particle!, root, the_q), particles)
    # foreach used instead of map because the results of the function calls are not needed.
    foreach(update_particle!, particles)
    nothing
end

# direct-sum step
function step!(particles::Vector{Particle})::Nothing
    for (p, q) ∈ combinations(particles, 2)
        # loop through pairwise combinations of particles,
        # equiallent to 2 for loops and an if statement skipping the case when both for loops give the same particle.

        fp::SVector{2, Float64} = calculate_force(p, q)
        # calculate the force on p with Newton's formula

        p.force_applied += fp
        q.force_applied += -fp # Newton's third law
    end
    update_particles!(particles)
    nothing
end
 
# ----- saving ------
# these functions are written generically to save and read a vector of any type

const delimiter::String = ";" # can't use comma because vectors use commas

function save!(file::IOStream, items::Vector{T}, data::Vector{Int64})::Nothing where {T}
    push!(data, length(items))
    # append the number of items to the data

    join(file, data, delimiter)
    # write the data to the file, joined by the delimiter

    write(file, "\n")
    # write a newline between data and the items

    fields::Tuple{Vararg{Symbol}} = fieldnames(Particle)
    for item ∈ items
        join(file, [getfield(item, f) for f ∈ fields], delimiter)
        # write each attribute of the item to the file, joined by delimiters
        # have to use comprehension instead of broadcast
        # because i can't garuntee that `length` will be defined for `item`

        write(file, "\n")
        # write a newline between each item
    end
    nothing
end

"""
    read(stream::IOStream, T::DataType)::Tuple{Tuple{Vararg{Int64}}, Vector{T}}

read data from IOStream as was written by `save!`
"""
function read(stream::IOStream, T::DataType)::Tuple{Tuple{Vararg{Int64}}, Vector}
    (data..., num_items) = parse.(Int64, tuple(string.(split(readline(stream), delimiter))...))
    # reading the data at the top of the file, the last number is always the number of items.
    # steps to read the numbers:
    # - read the first line of the file, where the numbers are
    # - split on delimiters
    # - convert the numbers to strings (split returns a `Vector` of `SubString`s, which `parse` can't read)
    # - convert the vector to a tuple
    # - parse each string in the tuple to an Int64
    # - unpack those into data and the number of items

    items::Vector{T} = Vector{T}(undef, num_items)
    # create an empty vector `num_items` long

    for (idx, line) ∈ enumerate(eachline(stream))
        item_data::Vector = [type(eval(Meta.parse(data))) for (type, data) ∈ zip(T.types, split(line, delimiter))]
        # `zip(T.types, split(line, delimiter))` matches up the attributes of `T`
        # to the attributes that were written to the file, split on the delimiter.
        # example output: `[(Int64, "22"), (String, "John"), (Float64, "89.5")]`

        # `Meta.parse` then parses the string to and expression,
        # which is then evaluated by `eval`
        # and then conveted to the type from T.types

        items[idx] = T(item_data...)
        # item_data then conveted to the type passed in by the caller and put in the vector.
    end
    return (data, items)
end

# ----- main function -----

function test_saving(particles::Vector{Particle})::Vector{Particle}
    fp = "data/save.txt"
    file = open(fp, "w")
    save!(file, particles, Vector{Int64}())
    close(file)

    file = open(fp, "r")
    data, ps = read(file, Particle)
    ps::Vector{Particle}
    close(file)
    return ps
end

function main()::Nothing
    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=4,
        edge_len=EDGE,
        mass_mean=10.0^24,
        mass_stddev=10.0^24,
        velocity_mean=10.0^24,
        velocity_stddev=0.0
        ))
    push!(particles, Particle(mass=10.0^30, fixed=true))

    # println(particles)
    # ps2 = test_saving(particles)
    # println()
    # println(ps2)
    # println(particles == ps2)

    t::Float64 = 0
    while !isempty(particles)
        start::DateTime = now()
        root = unsafe_from_just(BHTree(particles, SA[0.0, 0.0], 2 * EDGE))
        # construct the quadtree
        step!(particles, root)
        showparticles(particles)

        # contolling the timings:

        timetaken = convert(Millisecond, now() - start)
        # time taken to compute that frame, in milliseconds

        if (v = value(timetaken)) < Δt * 1000
            # if the timetaken is less that the target delta time
            sleeptime = Δt - v/1000
            # sleep for the differnce
            println("sleeping for $sleeptime")
            sleep(sleeptime)
            
        end
        t += Δt

    end
    nothing
end

main()
