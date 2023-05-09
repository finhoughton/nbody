using  LinearAlgebra: normalize, norm

struct Particle
    mass :: Float32
    pos :: Vector{Float32}
end

const G :: Float32 = 6.674e-11

function fast_inv_sqrt(x::Float32) :: Float32
    magic_num :: Int32 = 1597463007
    i = reinterpret(Int32, x)
    i = magic_num - (i >> 1) 
    y :: Float32 = reinterpret(Float32, i)
    return y * (1.5f0 - (x * 0.5f0 * y * y))
end

square(x::Number)::Number = x^2

newton_force_1(j::Particle, i::Particle)::Vector{Float32} = G * j.mass * i.mass * 1/(norm(j.pos-i.pos)^2) * normalize(j.pos - i.pos)

newton_force_2(j::Particle, i::Particle)::Vector{Float32} = G * j.mass * i.mass * (j.pos-i.pos) * 1/(norm(j.pos-i.pos)^ 3)

newton_force_3(j::Particle, i::Particle)::Vector{Float32} = G * j.mass * i.mass * (j.pos-i.pos) * (fast_inv_sqrt((j.pos[1] - i.pos[1])^2 + (j.pos[2] - i.pos[2])^2) ^ 3)

newton_force_4(j::Particle, i::Particle)::Vector{Float32} = G * j.mass * i.mass * (j.pos-i.pos) * (fast_inv_sqrt(sum(square.(j.pos-i.pos))) ^ 3)

function time_func(f :: Function, iterations :: Integer, args :: Vector{Particle}) :: Nothing
    for i in 1:iterations
        res = f(args...)
        if i == 1
            println("$f $res")
        end
    end
    return nothing
end

function main()
    a = Particle(100, [3.0f0, 4.0f0])
    b = Particle(120, [-6.0f0, 18.2f0])
    @time time_func(newton_force_1, 10000, [a, b])
    @time time_func(newton_force_2, 10000, [a, b])
    @time time_func(newton_force_3, 10000, [a, b])
    @time time_func(newton_force_4, 10000, [a, b])
    # just using normalize is faster than messing around with fast_inv_sqrt, thankfully.
end

main()
