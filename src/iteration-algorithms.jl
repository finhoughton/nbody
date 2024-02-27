using DataStructures

include("particle.jl")

const MAC::Float64 = 1 # smaller MAC = more accurate

struct BHTree
    particles::Vector{Particle}
    children::SVector{4, Maybe{BHTree}}
    NW::Maybe{BHTree}
    NE::Maybe{BHTree}
    SW::Maybe{BHTree}
    SE::Maybe{BHTree}
    centre::SVector{2, Float64}
    total_mass::Float64
    centre_of_mass::SVector{2, Float64}
    side_length::Float64

    function BHTree(
        particles::Vector{Particle},
        centre::SVector{2, Float64},
        side_length::Float64,
    )::Maybe{BHTree}
        
        len = length(particles)

        if side_length == 0
            error("0 length")
        end

        if len == 0
            return nothing
        end

        quadrants::Tuple{Vector{Particle}, Vector{Particle}, Vector{Particle}, Vector{Particle}} = ([], [], [], [])

        if len == 1
            total_mass = only(particles).mass
            centre_of_mass = only(particles).pos
            children = SVector{4, Maybe{BHTree}}([nothing, nothing, nothing, nothing])
        else

            # centre of mass = (m_1r_1 + m_2r_2 + ...)/(m_1 + m_2 + ...)
            total_mass::Float64 = 0
            total::SVector{2, Float64} = SA[0.0, 0.0]
            for particle ∈ particles
                total_mass += particle.mass
                total += particle.pos * particle.mass
                push_to_quadrant!(quadrants..., particle, centre)
            end
            centre_of_mass::SVector{2, Float64} = total / total_mass
            half_side_len::Float64 = 0.5 * side_length
            quarter_side_len::Float64 = 0.5 * half_side_len 
            centres = (
                centre + SVector(-quarter_side_len, quarter_side_len),  # NW
                centre + SVector(quarter_side_len, quarter_side_len),   # NE
                centre + SVector(-quarter_side_len, -quarter_side_len), # SW
                centre + SVector(quarter_side_len, -quarter_side_len))  # SE

                # recurive call creating the 4 children if there are more than 1 particle
                children = SVector{4, Maybe{BHTree}}([BHTree(quad, quad_centre, half_side_len) for (quad, quad_centre) ∈ zip(quadrants, centres)])
        end

        return Just(new(particles, children, children..., centre, total_mass, centre_of_mass, side_length))
    end
end

function Base.show(io::IO, x::BHTree)
    num_childs = length(filter((x -> is_something(x)), x.children))
    return print(io, "BHTree($(length(x.particles)) particles, centre $(x.centre), $(num_childs) children)")
end

Base.in(p::Particle, b::BHTree) = p in b.particles

function push_to_quadrant!(
    nws::Vector{Particle},
    nes::Vector{Particle},
    sws::Vector{Particle},
    ses::Vector{Particle},
    p::Particle,
    c::SVector{2, Float64}
    )::Nothing

    if p.pos[1] >= c[1]
        # it is west of the centre
        if p.pos[2] >= c[2]
            push!(nws, p)
            return nothing
        else
            push!(sws, p)
            return nothing
        end

    else
        # it is east of the centre
        if p.pos[2] >= c[2]
            push!(nes, p)
            return nothing
        else
            push!(ses, p)
            return nothing
        end
    end
    error("not pushed to vector")
end


function calculate_force(p::Particle, node::BHTree)::SVector{2, Float64}
    dinv = inv(EPS_SOFTENING + norm(p.pos - node.centre_of_mass) ^ 2)
    unit_v = normalize(node.centre_of_mass - p.pos) 
    F = unit_v * G * p.mass * (node.total_mass * dinv)
    println("force is $F")
    F
end
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
