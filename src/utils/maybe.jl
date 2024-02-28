
struct Just{T}
    __v::T
end

Maybe{T} = Union{Just{T}, Nothing}

function Base.:(==)(a::Just{T}, b::Just{T})::Bool where {T}
    a.__v == b.__v
end
 
"""
    from_maybe_with(default::A, func::(T -> A), value::Maybe{T})::A where {T, A}

If the `Maybe` value is `nothing`, the function returns the `default` value. 
Otherwise, it applies the function to the value inside the `Maybe` and returns the result.
"""
function from_maybe_with(default, func, value::Maybe{T}) where {T}
    is_nothing(value) ? default : func(value.__v)
end

"""
    from_maybe(default::T, value::Maybe{T})::T where {T}

If the Maybe is `nothing`, it returns the `default` value; otherwise, it returns the value contained in the `Maybe`.
"""
function from_maybe(default::T, value::Maybe{T})::T where {T}
    from_maybe_with(default, identity, value)
end

"""
    unsafe_from_just(value::Maybe{T})::T where {T}

Attempts to extract the value from a `Maybe`.
if it is `nothing`, an error is raised, otherwise the value is returned.
"""
function unsafe_from_just(value::Maybe{T})::T where {T}
    value.__v
end

"""
    branch_if_just(default::(A), func::(T -> A), value::Maybe{T})::A where {T, A}

if the value is `nothing`, call the defualt function and return the result.
Otherwise, extract the value from the `Just`, pass it to `func`, and return the result.
"""
function branch_if_just(default, func, value::Maybe{T}) where {T}
    is_nothing(value) ? default() : func(value.__v)
end

"""
    is_something(value::Maybe{T})::Bool where {T}

if a `Maybe` value contains a value, equivalent to `!is_nothing(value)`
"""
function is_something(value::Maybe{T})::Bool where {T}
    value ≢ nothing
end

"""
    is_nothing(value::Maybe{T})::Bool where {T}

if a `Maybe` value is `nothing`
"""
function is_nothing(value::Maybe{T})::Bool where {T}
    value ≡ nothing
end

"""
    liftM(func::(A, B, ... -> Z))::(Maybe{A}, Maybe{B}, ... -> Maybe{Z})

lift a function to be applicable to `Maybe` values. If any of the input values are `nothing`, `nothing` is returned.

# examples
```julia-repl
julia> maybe_add = liftM(+)
(::var"#lifted#"{typeof(+)}) (generic function with 1 method)
julia> a = maybe_add(Just(10), Just(4), Just(100))
Just{Int64}(114)

julia> nothing === maybe_add(Just(10), Just(4), Just(100), nothing)
true
```

"""
function liftM(func)
    function lifted(vs...)
        any((is_nothing(v) for v ∈ vs)) ? nothing : Just(func((v.__v for v ∈ vs)...))
    end
end


"""
    drop_nothings(xs::Vector{Maybe{T}})::Vector{T} where {T}

remove nothings from a vector, keeping only the values inside `Maybe`s
"""
function drop_nothings(xs::AbstractVector{Maybe{T}})::Vector{T} where {T}
    [x.__v for x ∈ xs if is_something(x)]
end

nothing
