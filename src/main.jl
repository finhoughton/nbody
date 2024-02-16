using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Combinatorics: combinations
using Random
using GLMakie
using PartialFunctions
using Base.Threads

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
    update_particle!.(particles)
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

function get_unnamed_filename()::Int
    files::Vector{String} = filter(n -> endswith(n, ".txt") && startswith(n, "untitled"), readdir("data"))
    # get a list of files in the data/ directory, and filter by txt files starting with "untitled"
    if isempty(files)
        return 1
    end
    filenums::Vector{Int} = map((parse $ Int) ∘ String ∘ collect ∘ (Iterators.takewhile $ isdigit) ∘ (Iterators.dropwhile $ !isdigit), files)
    # extract the numbers from the untitled files
    return max(filenums...) + 1
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

function step_gui!(xs, ys, particles::Vector{Particle})::Nothing
    # root = BHTree(particles, SA[0.0, 0.0], 2 * EDGE) |> unsafe_from_just
    # step!(particles, root)
    start = now()
    step!(particles)
    timetaken = convert(Millisecond, now() - start)
    xs[] = [p.pos[1] for p ∈ particles]
    ys[] = [p.pos[2] for p ∈ particles]
    if (v = value(timetaken)) < Δt * 1000
        # if the timetaken is less that the target delta time
        sleeptime = Δt - v/1000
        # sleep for the differnce
        # println("sleeping for $sleeptime")
        sleep(sleeptime)
        
    end
end

function main()::Nothing

    # adding some particles
    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=8,
        edge_len=EDGE,
        mass_mean=10.0^24,
        mass_stddev=10.0^23,
        velocity_mean=0.0,
        velocity_stddev=5 * 10.0^8,
        ))
    push!(particles, Particle(mass=10.0^30, fixed=true))

    for (i, p) ∈ enumerate(particles)
        p.id = i
    end

    t::Float64 = 0
    iteration_number::Int64 = 0

    # -- GUI setup -- 

    fig = Figure()
    display(fig)

    ax = Axis(fig[1, 1])

    xs::Observable{Vector{Float64}} = Observable([p.pos[1] for p ∈ particles])
    ys::Observable{Vector{Float64}} = Observable([p.pos[2] for p ∈ particles])
    scatter!(ax, xs, ys; color=:blue, marker=:circle, markersize=10)
    ax.title = "Unititled Simulation"
    # xlims!(ax, -EDGE, EDGE)
    # ylims!(ax, -EDGE, EDGE)

    fig[2, 1] = buttons = GridLayout(tellwidth = false)

    start_stop_button = Button(buttons[1, 1]; label = "start/stop", tellwidth = false)
    
    is_running = Observable(false)
    
    function update_callback(t)
        iteration_number += 1
        step_gui!(xs, ys, particles)
    end
    
    on(start_stop_button.clicks) do clicks
        is_running[] = !is_running[] # switch the state of is_running
        
        if is_running[]
            global timer
            timer = Timer(update_callback, 0, interval=Δt)
        else
            close(timer)
        end
    end

    save_button = Button(buttons[1, 2]; label = "Save simulaton")
    savename_tb = Textbox(buttons[1, 3], placeholder = "Enter save name", tellwidth = false)
    load_button = Button(buttons[1, 4]; label = "Load simulaton")
    stored_filename = ""

    on(savename_tb.stored_string) do s
        stored_filename = string(s)
    end


    on(save_button.clicks) do clicks
        if is_running[]
            close(timer)
        end

        filename = stored_filename == "" ? string("untitled", get_unnamed_filename()) : stored_filename

        f = open(string("data/", filename, ".txt"), "w")
        save!(f, particles, Vector{Int64}())
        close(f)
        savename_tb.stored_string = ""
    end


    while true
        # println("$iteration_number, $(is_running[]), $xs, $ys")
        sleep(0.01)
    end
    nothing
end

main()
