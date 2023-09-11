
struct Maybe{A}
    value::Union{A, Nothing}
end

function to_maybe(f, default::A, value::A)::Maybe{A}
    f(value) ? Maybe(value) : default
end

function maybe(f, default::B, o::Maybe{A})::B
    is_nothing(o) ? default : f(o.value)
end

function from_maybe( default::A, o::Maybe{A})::A
    maybe(identity, default, o)
end

function unsafe_from_maybe(o::Maybe{A})::A
    is_nothing(o) ? error("unsafe_from_maybe recieved nothing value") : o.value
end

function has_value(o::Maybe{A})::Bool
    o.value !== nothing
end

function is_nothing(o::Maybe(A))::Bool
    o.value === nothing
end

