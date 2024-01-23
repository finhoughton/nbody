using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Combinatorics: combinations
using Random
using GLMakie
using PartialFunctions

include("iteration-algorithms.jl")

# order of include():
# main
# iteration algorithms
# particle
# utils

# ------ displaying ------
# using Plots
# gr()
# function showparticles(particles::Vector{Particle})::Nothing
#     positions::Vector{Vector{Float64}} = [p.pos for p in particles]
#     xs::Vector{Float64} = getindex.(positions, 1)
#     ys::Vector{Float64} = getindex.(positions, 2)
#     p::Plots.Plot = scatter(xs, ys, legend=false, showaxis=false, background_colour=:black)
#     xlims!(p, X_LIMITS...)
#     ylims!(p, Y_LIMITS...)
#     display(p)
#     nothing
# end

# ------ iteration ------



# barnes-hut
function step_particle!(root::BHTree, the_q::Queue{BHTree}, p::Particle)::Nothing 
    function my_enqueue!(x)
        println("enqueuing $x")
        enqueue!(the_q, x)
    end
    my_enqueue!(root)
    # enqueue the root node
    while !isempty(the_q)
        println("the_q $the_q")
        current::BHTree = dequeue!(the_q)
        # dequeue a node
        distance_to_centre::Float64 = norm(p.pos - current.centre)
        # calculate the distance to the centre to be used in MAC calculations
        children = drop_nothings(current.children)
        if current.side_length < MAC * distance_to_centre || isempty(children)
            # the ratio is not greater than the MAC, use this quadrant to appriximate the force on the particle
            if length(current.particles) != 1 || only(current.particles) != p
                p.force_applied += calculate_force(p, current)
            end
        else
            # the ratio is greater than the MAC, enqueue the current node's children
            my_enqueue!.(children)
        end
    end
    nothing
end

# barnes-hut step
function step!(particles::Vector{Particle}, root::BHTree)::Nothing
    the_q = Queue{BHTree}()
    # create the queue used in bh algorithm
    foreach(step_particle! $ (root, the_q), particles)
    # foreach used instead of map because the results of the function calls are not needed.
    f = [p.force_applied for p in particles]
    println("before update_particle!, $f")
    foreach(update_particle!, particles)
    f = [p.force_applied for p in particles]
    println("after update_particle!, $f")
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

split_delim(l) = split(l, delimiter)

function save!(file::IOStream, items::Vector{T}, data::Vector{Int64})::Nothing where {T}
    push!(data, length(items))
    # append the number of items to the data

    join(file, data, delimiter)
    # write the data to the file, joined by the delimiter

    write(file, "\n")
    # write a newline between data and the items

    fields::Tuple{Vararg{Symbol}} = fieldnames(Particle)
    for item ∈ items
        join(file, map(getfield $ item, fields), delimiter)
        # write each attribute of the item to the file, joined by delimiters
        # have to use map instead of broadcast
        # because `length` may not be defined for `item`

        write(file, "\n")
        # write a newline between each item
    end
    nothing
end

function read!(stream::IOStream, T::DataType)::Tuple{Tuple{Vararg{Int64}}, Vector}
    (data..., num_items) = stream |> readline |> split_delim |> (x -> parse.(Int64, x)) |> Tuple
    # reading the data at the top of the file, the last number is always the number of items.
    # steps to read the numbers:
    # - read the first line of the file, where the numbers are, to a string
    # - split the string on delimiters
    # - parse each substring  to an Int64
    # - convert to a tuple
    # - unpack those into data and the number of items

    items::Vector{T} = Vector{T}(undef, num_items)
    # create an empty vector `num_items` long

    for (idx, line) ∈ enumerate(eachline(stream))
        item_data::Vector = [(data |> Meta.parse |> eval |> type) for (type, data) ∈ (line |> split_delim |> zip $ T.types)]
        # `line |> split_delim |> zip $ T.types` splits the line on delimiters
        # then matches up the attributes of `T` to the attributes that were written to the file
        # example output: `[(Int64, "22"), (String, "John"), (Float64, "89.5")]`

        # `Meta.parse` then parses the string to an expression,
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

function step_gui!(positions::Observable{Vector{SVector{2, Float64}}}, particles::Vector{Particle})::Nothing
    positions[] = [p.pos for p ∈ particles]
end

function main()::Nothing

    # adding some particles
    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=5,
        edge_len=EDGE,
        mass_mean=10.0^24,
        mass_stddev=10.0^24,
        velocity_mean=10.0^24,
        velocity_stddev=0.0
        ))
    push!(particles, Particle(mass=10.0^30, fixed=true))

    for (i, p) ∈ enumerate(particles)
        p.id = i
    end

    t::Float64 = 0

    # -- GUI setup -- 

    fig = Figure()
    display(fig)

    ax = Axis(fig[1, 1])

    positions::Observable{Vector{SVector{2, Float64}}} = Observable([p.pos for p ∈ particles])
    scatter!(ax, positions; color=:blue, marker=:circle, markersize=2)
    ax.title = "Unititled Simulation"
    xlims!(ax, -EDGE, EDGE)
    ylims!(ax, -EDGE, EDGE)

    # -- start/stop button --

    start_stop = Button(fig[2,1]; label = "start/stop", tellwidth = false)

    is_running = Observable(false)

    on(start_stop.clicks) do clicks
        is_running[] = !is_running[] # switch to state of is_running
    end

    on(start_stop.clicks) do clicks
        @async while is_running[]
            isopen(fig.scene) || break # stop calculations if closed window
            step_gui!(positions, particles)
            yield()
        end
    end

    # -- main loop -- 
    i = 0
    while true
        start::DateTime = now()
        println("i $i")
        for p ∈ particles
            println(p.pos)
        end
        i += 1
        root = BHTree(particles, SA[0.0, 0.0], 2 * EDGE) |> unsafe_from_just
        # construct the quadtree
        f = [p.force_applied for p in particles]
        println("just after bhtree, $f")
        step!(particles, root)
        f = [p.force_applied for p in particles]
        println("just after step, $f")
        # showparticles(particles)



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
