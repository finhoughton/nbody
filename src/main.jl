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
    xs.val = [p.pos[1] for p ∈ particles]
    ys.val = [p.pos[2] for p ∈ particles]
    notify(xs)
    notify(ys)
    nothing
end

function main()::Nothing

    # -- particle setup --

    particles::Vector{Particle} = []
    append!(particles, random_particles(
        n=50,
        edge_len=EDGE,
        mass_mean=10.0^25,
        mass_stddev=5 * 10.0^24,
        velocity_stddev=5 * 10.0^8,
    ))
    push!(particles, Particle(mass=10.0^30, fixed=false))

    for (i, p) ∈ enumerate(particles)
        p.id = i
    end

    iteration_num::Int64 = 0

    # -- GUI setup -- 

    fig = Figure()
    display(fig)
    resize!(fig.scene, 1000, 1000)

    ax = Axis(fig[1, 2])

    xs::Observable{Vector{Float64}} = Observable([p.pos[1] for p ∈ particles])
    ys::Observable{Vector{Float64}} = Observable([p.pos[2] for p ∈ particles])
    scatter!(ax, xs, ys; color=:blue, marker=:circle, markersize=10)
    ax.title = "Unititled Simulation"

    fig[2, 2] = menu = GridLayout(tellwidth=false)

    # -- stats -- 

    num_particles_lb = Label(menu[1, 1], "particles: 0")
    ups_lb = Label(menu[2, 1], "UPS: 0")

    prev_iteration_num::Int64 = 0

    function update_stats!(t::Timer)::Nothing
        Δi = iteration_num - prev_iteration_num
        prev_iteration_num = iteration_num
        ups_lb.text = "UPS: $Δi"
        num_particles_lb.text = "Particles: $(length(particles))"
        nothing
    end
    Timer(update_stats!, 0, interval=1)
    # update the stats every second

    # -- barnes-hut toogle --

    Label(menu[1, 2], "Use Barnes-Hut?")
    bh_toogle = Toggle(menu[1, 3]; active=false)

    use_bh::Bool = false
    on(bh_toogle.active) do a
        string("Switched to ", a ? "Barnes-Hut" : "direct sum") |> println
        use_bh = a
    end

    # -- MAC box --

    Label(menu[2, 2], "MAC: ")
    mac_tb = Textbox(menu[2, 3], placeholder="1", validator=Float64)

    on(mac_tb.stored_string) do s
        mac = parse(Float64, s)
        println("mac set to $mac")
    end

    # -- particle max speed --

    speed_lb = Label(menu[5, 1], "particle max speed is:\nInf")
    max_speed_tb = Textbox(menu[5, 2], placeholder="Inf", validator=(s -> (x = tryparse(Float64, s)) ≢ nothing && x > 0))

    particle_max_speed::Float64 = Inf
    on(max_speed_tb.stored_string) do s
        particle_max_speed = parse(Float64, s)
        println("particle max speed set to $particle_max_speed")
        speed_lb.text = "particle max speed is:\n$particle_max_speed"
    end

    # -- speed slider, controlls actually stepping the simulation --

    function step_sim!(t::Timer)
        if use_bh
            root::BHTree = make_bh(particles)
            step!(root, particles)
        else
            step!(particles)
        end
        foreach(update_particle!$particle_max_speed, particles)
        iteration_num += 1
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
        println("slider set to $v")
    end

    function pause_sim!()
        set_close_to!(speed_sl, 0)
    end

    # -- save button and text box --

    save_button = Button(menu[1, 4]; label="Save simulaton to: untited.txt")
    savename_tb = Textbox(menu[2, 4], placeholder="unitited", tellwidth=false)

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

    load_button = Button(menu[1, 5]; label="Load simulaton from: ")
    loadname_tb = Textbox(menu[2, 5], placeholder="untited", tellwidth=false)

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

    generate_button = Button(menu[4, 1]; label="Generate particless")

    Label(menu[3, 2], "n:")
    num_particles_tb = Textbox(menu[4, 2], placeholder="0", validator=Int64)

    num_particles::Int64 = 0
    on(num_particles_tb.stored_string) do s
        num_particles = parse(Int64, s)
    end

    Label(menu[3, 3], "mass mean:")
    mass_mean_tb = Textbox(menu[4, 3], placeholder="1e24", validator=Float64)
    mass_mean::Float64 = 1e24
    on(mass_mean_tb.stored_string) do s
        mass_mean = parse(Float64, s)
    end

    Label(menu[3, 4], "mass stdddev:")
    mass_stddev_tb = Textbox(menu[4, 4], placeholder="0", validator=Float64)
    mass_stddev::Float64 = 0
    on(mass_stddev_tb.stored_string) do s
        mass_stddev = parse(Float64, s)
    end

    Label(menu[3, 5], "velocity stdddev:")
    velocity_stddev_tb = Textbox(menu[4, 5], placeholder="0", validator=Float64)
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
        sleep(Δt/100)
    end
    close.(timers)
    nothing
end

main()
