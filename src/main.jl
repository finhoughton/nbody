using LinearAlgebra
using StaticArrays
using Combinatorics: combinations
using Random
using PartialFunctions
using GLMakie

include("iteration-algorithms.jl")
include("saving.jl")

# order of include():
# main
# - saving
# - iteration algorithms
# - - particle
# - - - utils

function update_gui!(xs, ys, particles::Vector{Particle})::Nothing
    xs[] = [p.pos[1] for p ∈ particles]
    ys[] = [p.pos[2] for p ∈ particles]
    nothing
end

function main()::Nothing

    # -- particle setup --

    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=300,
        edge_len=EDGE,
        mass_mean=10.0^25,
        mass_stddev=5 * 10.0^24,
        velocity_stddev=5 * 10.0^8,
    ))
    push!(particles, Particle(mass=10.0^30, fixed=false))

    for (i, p) ∈ enumerate(particles)
        p.id = i
    end

    t::Float64 = 0

    # -- GUI setup -- 

    fig = Figure()
    display(fig)
    resize!(fig.scene, 1000, 500)

    ax = Axis(fig[1, 2])

    xs::Observable{Vector{Float64}} = Observable([p.pos[1] for p ∈ particles])
    ys::Observable{Vector{Float64}} = Observable([p.pos[2] for p ∈ particles])
    scatter!(ax, xs, ys; color=:blue, marker=:circle, markersize=10)
    ax.title = "Unititled Simulation"

    fig[2, 2] = buttons = GridLayout(tellwidth=false)

    # -- barnes-hut toogle --

    bh_toogle = Toggle(buttons[1, 2]; active=false)
    bh_label = Label(buttons[1, 1], "Use Barnes-Hut?")

    use_bh::Bool = false
    on(bh_toogle.active) do a
        string("Switched to ", a ? "Barnes-Hut" : "direct sum") |> println
        use_bh = a
    end

    # -- speed slider, controlls actually stepping the simulation --

    function step_sim!(t::Timer)
        if use_bh
            root::BHTree = make_bh(particles)
            step!(root, particles)
        else
            step!(particles)
        end
    end

    # each timer iterates the simulation every Δt, more timers = faster
    timers::Vector{Timer} = []

    function set_num_timers!(n::Int)::Nothing
        l = length(timers) # the current number of timers
        Δn = n - l
        if Δn == 0
            return nothing
        elseif Δn > 0
            # we need to create more timers
            ΔΔt = Δt / Δn
            [Timer((_::Timer) -> push!(timers, Timer(step_sim!, 0, interval=Δt)), (i - 1) * ΔΔt) for i ∈ 1:Δn]
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

    # -- save button and text box --

    save_button = Button(buttons[1, 3]; label="Save simulaton to: untited.txt")
    savename_tb = Textbox(buttons[2, 3], placeholder="unitited", tellwidth=false)

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
        save!(f, particles, Int64[])
        close(f)
    end

    # -- load button and text box --

    load_button = Button(buttons[1, 4]; label="Load simulaton from: ")
    loadname_tb = Textbox(buttons[2, 4], placeholder="untited", tellwidth=false)

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
            Timer(function (t::Timer) # reset the load button in 2 seconds
                    load_button.label = "Load simulation from: "
                    load_button.labelcolor = :black
                end, 2)
        end
    end

    # -- buttons for generating particles

    generate_button = Button(buttons[4, 1]; label="Generate particless")

    Label(buttons[3, 2], "n:")
    num_particles_tb = Textbox(buttons[4, 2], placeholder="0", validator=Int64)

    num_particles::Int64 = 0
    on(num_particles_tb.stored_string) do s
        num_particles = parse(Int64, s)
    end

    Label(buttons[3, 3], "mass mean:")
    mass_mean_tb = Textbox(buttons[4, 3], placeholder="1e24", validator=Float64)
    mass_mean::Float64 = 1e24
    on(mass_mean_tb.stored_string) do s
        mass_mean = parse(Float64, s)
    end

    Label(buttons[3, 4], "mass stdddev:")
    mass_stddev_tb = Textbox(buttons[4, 4], placeholder="0", validator=Float64)
    mass_stddev::Float64 = 0
    on(mass_stddev_tb.stored_string) do s
        mass_stddev = parse(Float64, s)
    end

    Label(buttons[3, 5], "velocity stdddev:")
    velocity_stddev_tb = Textbox(buttons[4, 5], placeholder="0", validator=Float64)
    velocity_stddev::Float64 = 0
    on(velocity_stddev_tb.stored_string) do s
        velocity_stddev = parse(Float64, s)
    end

    on(generate_button.clicks) do _
        ps::Vector{Particle} = random_particles(
            n=num_particles,
            edge_len=EDGE,
            mass_mean=mass_mean,
            mass_stddev=mass_stddev,
            velocity_stddev=velocity_stddev
        )
        push!(particles, ps...)
    end

    # -- main loop --

    while isopen(fig.scene)
        update_gui!(xs, ys, particles)
        sleep(Δt)
    end
    nothing
end

main()
