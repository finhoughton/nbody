using LinearAlgebra
using StaticArrays
using Dates: now, value, DateTime, Millisecond
using Combinatorics: combinations
using Random
using PartialFunctions
using GLMakie

include("iteration-algorithms.jl")
include("saving.jl")

# order of include():
# main
# iteration algorithms
# particle
# utils

function update_gui!(xs, ys, particles::Vector{Particle})::Nothing
    xs[] = [p.pos[1] for p ∈ particles]
    ys[] = [p.pos[2] for p ∈ particles]
    nothing
end

function main()::Nothing

    # adding some particles
    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=30,
        edge_len=EDGE,
        mass_mean=10.0^25,
        mass_stddev=5*10.0^24,
        velocity_stddev=5*10.0^8,
        ))
    push!(particles, Particle(mass=10.0^29, fixed=false))

    for (i, p) ∈ enumerate(particles)
        p.id = i
    end

    t::Float64 = 0

    # -- GUI setup -- 

    fig = Figure()
    display(fig)

    ax = Axis(fig[1, 2])

    xs::Observable{Vector{Float64}} = Observable([p.pos[1] for p ∈ particles])
    ys::Observable{Vector{Float64}} = Observable([p.pos[2] for p ∈ particles])
    scatter!(ax, xs, ys; color=:blue, marker=:circle, markersize=10)
    ax.title = "Unititled Simulation"

    fig[2, 2] = buttons = GridLayout(tellwidth = false)
    
    is_running = Observable(false)
    
    function update_callback(t::Timer)
        step_sim!(particles)
        update_gui!(xs, ys, particles)
    end


    # each timer iterates the simulation every Δt, more timers = faster
    timers::Vector{Timer} = []

    function set_num_timers!(n::Int)::Nothing
        l = length(timers)
        Δn = n - l
        if Δn == 0
            return nothing
        elseif Δn > 0
            # we need to create  more timers
            ΔΔt = Δt / Δn
            [Timer((_::Timer) -> push!(timers, Timer(update_callback, 0, interval=Δt)), (i - 1) * ΔΔt) for i ∈ 1:Δn]
        else
            # we need to close some timers
            [close(pop!(timers)) for _ ∈ 1:abs(Δn)]
        end
        return nothing
    end

    speed_sl = Slider(fig[1, 1], range=0:1:20, startvalue=0, horizontal=false)
 
    on(speed_sl.value) do v
        set_num_timers!(v)
    end

    function pause_sim!()
        set_close_to!(speed_sl, 0)
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
        pause_sim!()

        filename = (save_file == "") ? string("untitled", get_untitled_filename_number()) : save_file

        ax.title = filename

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
        if isfile(filename)
            pause_sim!()
            ax.title = load_file
            f = open(filename, "r")
            _, ps = read!(f, Particle)
            ps::Vector{Particle}
            particles = ps
            update_gui!(xs, ys, particles)
        else
            load_button.label = "file does not exist!"
            load_button.labelcolor = :red
            Timer(function(t::Timer) # reset the load button in 2 seconds
                load_button.label = "Load simulation from: "
                load_button.labelcolor = :black
            end, 2)
        end
    end



    while isopen(fig.scene)
        # println("$iteration_number, $(is_running[]), $xs, $ys")
        sleep(Δt)
    end
    nothing
end

main()
