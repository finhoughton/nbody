using DataStructures
using StaticArrays

include("particle.jl")

mac::Float64 = 1 # smaller MAC = more accurate

mutable struct BHTree
    particles::Vector{Particle}
    children::Vector{Maybe{BHTree}} # NW, NE, SW, SE
    centre::SVector{2,Float64}
    total_mass::Float64
    centre_of_mass::SVector{2,Float64}
    side_length::Float64
end

empty_bh(edge::Float64, centre::SVector{2,Float64}) = BHTree(Vector{Particle}(), Vector{Maybe{BHTree}}([nothing, nothing, nothing, nothing]), centre, 0.0, SA[0.0, 0.0], edge)

singleon_bh(edge::Float64, centre::SVector{2,Float64}, p::Particle) = BHTree([p], Vector{Maybe{BHTree}}([nothing, nothing, nothing, nothing]), centre, 0.0, SA[0.0, 0.0], edge)

function Base.push!(bh::BHTree, p::Particle)::Maybe{BHTree}
    push!(bh.particles, p)
    quad = 1 + Int(p.pos[1] > bh.centre[1]) + 2 * Int(p.pos[2] < bh.centre[2])
    if !is_nothing(bh.children[quad])
        push!(unsafe_from_just(bh.children[quad]), p)
        return nothing
    end

    if quad == 1
        direction = SA[-1, 1]
    elseif quad == 2
        direction = SA[1, 1]
    elseif quad == 3
        direction = SA[-1, -1]
    elseif quad == 4
        direction = SA[1, -1]
    else
        error("invalid quadrant")
    end
    new_centre::SVector{2, Float64} = bh.centre + 0.25 * (direction .* [bh.side_length, bh.side_length])
    child::Maybe{BHTree} = Just(singleon_bh(bh.side_length / 2, new_centre, p))
    bh.children[quad] = child
    Just(bh)
end

function calculate_centres_of_mass!(bh::BHTree)::Nothing
    total_mass::Float64 = sum(p.mass for p ∈ bh.particles) 
    total_mass_pos::SVector{2, Float64} = sum(p.mass * p.pos for p ∈ bh.particles)
    bh.total_mass = total_mass
    bh.centre_of_mass = total_mass_pos / total_mass
    foreach(calculate_centres_of_mass!, drop_nothings(bh.children))
    # calculate and store the centeres of mass of each of this node's child nodes
    nothing
end

function calculate_force(p::Particle, node::BHTree)::SVector{2,Float64}
    normalize(node.centre_of_mass - p.pos) * G * p.mass * (node.total_mass * inv(EPS_SOFTENING + norm(p.pos - node.centre_of_mass)^2))
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
        if current.side_length < mac * distance_to_centre || all(is_nothing, current.children)
            # the ratio is not greater than the MAC, use this quadrant to appriximate the force on the particle
            if length(current.particles) != 1 || only(current.particles) != p
                p.force_applied += calculate_force(p, current)
            end
        else
            # the ratio is greater than the MAC, enqueue the current node's children
            foreach(enqueue!$the_q, drop_nothings(current.children))
        end
    end
    nothing
end

# barnes-hut step
function step!(root::BHTree, particles::Vector{Particle})::Nothing
    the_q::Queue{BHTree} = Queue{BHTree}()
    # create the queue used in bh algorithm
    foreach(step_particle!$(root, the_q), particles)
    foreach(update_particle!, particles)
    # foreach used instead of map because the results of the function calls are not needed.
    nothing
end

# direct-sum step
function step!(particles::Vector{Particle})::Nothing
    for (p, q) ∈ combinations(particles, 2)
        # loop through pairwise combinations of particles,
        # equiallent to 2 for loops and an if statement skipping the case when both for loops give the same particle.

        fp::SVector{2,Float64} = calculate_force(p, q)
        # calculate the force on p with Newton's formula

        p.force_applied += fp
        q.force_applied += -fp # Newton's third law
    end
    update_particle!.(particles)
    nothing
end

function make_bh(ps::Vector{Particle})::BHTree
    e = 2.1 * max(maximum(abs, (p.pos[1] for p ∈ ps)), maximum(abs, (p.pos[2] for p ∈ ps)))
    # make sure every particle is within the root node
    bh = empty_bh(e, SA[0.0, 0.0])
    foreach(push!$bh, ps)
    calculate_centres_of_mass!(bh)
    return bh
end
