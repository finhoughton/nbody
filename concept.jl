using Combinatorics: combinations, Combinations
using LinearAlgebra: normalize, norm
using Dates: now, value, DateTime, Millisecond
using Plots
gr()

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 30

const EPS_SOFTENING::Float64 = 1e7
# stop forcess becoming too big when objects are very close

let id::Int = 0
    mutable struct Particle
        id::Int
        mass::Float64
        pos::Vector{Float64}
        v::Vector{Float64}
        force_applied::Vector{Float64}
        fixed::Bool
        function Particle(
            mass::Float64,
            pos::Vector{Float64},
            v::Vector{Float64}=zeros(Float64, 2),
            fixed::Bool=false
        )::Particle
            if size(pos, 1) != 2 && size(v, 1) != 2
                error("Particle only supports 2D positions and velocity")
            elseif mass <= 0
                error("particle mass must be positive")
            elseif fixed && norm(v) != 0
                error("fixed is incomaptable with velocity.")
            end
            x += 1
            new(x, mass, pos, v, zeros(Float64, 2), fixed)
        end
    end
end

function apply_force!(force::Vector{Float64}, particle::Particle)::Nothing
    particle.force_applied += force
    nothing
end

function force_pairwise!(j::Particle, i::Particle)::Nothing
    distance_inv::Float64 = 1/(norm(j.pos-i.pos) + EPS_SOFTENING)
    f_ji::Vector{Float64} = normalize(j.pos - i.pos) * -G * j.mass * i.mass * distance_inv * distance_inv
    f_ij::Vector{Float64} = -f_ji
    apply_force!(f_ji, j)
    apply_force!(f_ij, i)
    nothing
end

function step!(particle::Particle, dt::Float64) :: Nothing
    if particle.fixed
        # don't compute any force for fixed particles.
        return nothing
    end
    f::Vector{Float64} = particle.force_applied
    particle.force_applied = zeros(Float64, 2)
    da::Vector{Float64} = f/particle.mass
    dv::Vector{Float64} = da * dt
    particle.v += dv
    ds::Vector{Float64} = particle.v * dt
    particle.pos += ds
    nothing
end

function step!(particles::Vector{Particle}, dt::Float64)::Nothing
    step!.(particles, dt)
    nothing
end

function velocity_arrow!(particle::Particle)::Nothing
    # @assert particle.pos[1]/particle.pos[2] ≈ pos_norm[1]/pos_norm[2]
    arrow_start::Vector{Float64} = particle.pos
    arrow_end::Vector{Float64} = arrow_start + normalize(particle.v) * DIST/20
    plot!(arrow_start, arrow_end, arrow=true, color=:white, linewidth=2, label="")  
    nothing
end

function showparticles(particles::Vector{Particle})::Nothing

    positions::Vector{Vector{Float64}} = [p.pos for p in particles]

    xs::Vector{Float64} = getindex.(positions, 1)
    ys::Vector{Float64} = getindex.(positions, 2)
    p = scatter(xs, ys, legend=false, reuse=true, show=true, background_colour=:black)
    # velocity_arrow!.(particles)
    xlims!(-10 * DIST, 10 * DIST)
    ylims!(-10 * DIST, 10 * DIST)
    display(p)
    nothing
end

function random_particle() :: Particle
    m::Float64 = M_EARTH * (rand() + 0.5)
    pos::Vector{Float64} = (rand(Float64, 2) * 2 - ones(2)) * DIST * 4
    v::Vector{Float64} = (rand(Float64, 2) * 2 - ones(2)) * DIST * 0.5
    Particle(m, pos, v)
end

function main()
    N::Int = 12
    particles::Vector{Particle} = [random_particle() for i in 1:N]
    m_1::Particle=Particle(M_EARTH*50, [-5 * DIST, 0], zeros(Float64, 2), true)
    m_2::Particle=Particle(M_EARTH*50, [5 * DIST, 0], zeros(Float64, 2), true)
    push!(particles, m_1)
    push!(particles, m_2)
    t::Float64 = 0
    fps::Int = 30
    dt::Float64 = 1/fps

    while true
        start::DateTime = now()
        for pair in combinations(particles, 2)
            force_pairwise!(pair...)
        end

        step!(particles, dt)
        showparticles(particles)
        # for particle in particles
        #     println(round.(particle.v, digits=4))
        # end

        timetaken = convert(Millisecond, now() - start)
        if (v = value(timetaken)) < dt * 1000
            a = dt - v/1000
            println("sleeping for $a")
            sleep(a)
        end
        t += dt
    end
end

main()