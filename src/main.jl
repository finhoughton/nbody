using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Combinatorics: combinations
using Random
using GLMakie
using PartialFunctions
using Base.Threads

include("iteration-algorithms.jl")
include("saving.jl")

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

function step_sim!(particles::Vector{Particle})::Nothing
    # root = BHTree(particles, SA[0.0, 0.0], 2 * EDGE) |> unsafe_from_just
    # step!(particles, root)
    start = now()
    step!(particles)
    timetaken = convert(Millisecond, now() - start)
    if (v = value(timetaken)) < Δt * 1000
        # if the timetaken is less that the target delta time
        sleeptime = Δt - v/1000
        # sleep for the differnce
        # println("sleeping for $sleeptime")
        sleep(sleeptime)        
    end
end

function update_gui!(xs, ys, particles::Vector{Particle})::Nothing
    xs[] = [p.pos[1] for p ∈ particles]
    ys[] = [p.pos[2] for p ∈ particles]
    nothing
end

function main()::Nothing

    # adding some particles
    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=15,
        edge_len=EDGE,
        mass_mean=10.0^24,
        mass_stddev=10.0^23,
        velocity_mean=0.0,
        velocity_stddev=5 * 10.0^8,
        ))
    push!(particles, Particle(mass=10.0^31, fixed=true))

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

    fig[2, 1] = buttons = GridLayout(tellwidth = false)

    start_stop_button = Button(buttons[1, 1]; label = "start", tellwidth = false)
    
    is_running = Observable(false)
    
    function update_callback(t)::Nothing
        iteration_number += 1
        step_sim!(particles)
        update_gui!(xs, ys, particles)
        nothing
    end
    
    function start_sim!()::Nothing
        start_stop_button.label = "Stop"
        is_running[] = true
        global timer = Timer(update_callback, 0, interval=Δt)
        nothing
    end

    function stop_sim!()::Nothing
        close(timer)
        is_running[] = false
        start_stop_button.label = "Start"
        nothing
    end

    on(start_stop_button.clicks) do _
        println("start stop clicked, is running: $(is_running[])")
        if is_running[]
            stop_sim!()
        else
            start_sim!()
        end
    end

    save_button = Button(buttons[1, 2]; label = "Save simulaton to: untited.txt")
    savename_tb = Textbox(buttons[2, 2], placeholder = "unitited", tellwidth = false)
    load_button = Button(buttons[1, 3]; label = "Load simulaton from: ")
    loadname_tb = Textbox(buttons[2, 3], placeholder = "untited", tellwidth = false)


    # -- save button and text box --
    save_file = ""
    on(savename_tb.stored_string) do s
        save_file = string(s)
        save_button.label = string("Save simulation to: ", s, ".txt")
    end


    on(save_button.clicks) do _
        if is_running[]
            stop_sim!()
        end

        filename = (save_file == "") ? string("untitled", get_untitled_filename_number()) : save_file

        f = open(string("data/", filename, ".txt"), "w")
        save!(f, particles, Vector{Int64}())
        close(f)
    end

    # -- load button and text box --
    load_file = ""
    on(loadname_tb.stored_string) do s
        load_file = string(s)
        load_button.label = string("Load simulation from: ", s, ".txt")
    end

    on(load_button.clicks) do _
        filename = string("data/", load_file, ".txt")
        if !isfile(filename)
            load_button.label = "file does not exist!"
            for _ in 1:8
                load_button.labelcolor = :red
                sleep(0.3)
                load_button.labelcolor = :black
            end
            load_button.label = "Load simulation from: "
        else
            stop_sim!()
            f = open(filename, "r")
            _, ps = read!(f, Particle)
            ps::Vector{Particle}
            particles = ps
            update_gui!(xs, ys, particles)
        end
    end



    while true
        # println("$iteration_number, $(is_running[]), $xs, $ys")
        sleep(Δt)
    end
    nothing
end

main()
