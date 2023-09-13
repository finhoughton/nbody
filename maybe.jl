
struct Just{T}
    _v::T
end

Maybe{T} = Union{Just{T}, Nothing}

function ==(a::Just{T}, b::Just{T})::Bool where {T}
    a._v == b._v
end

function ==(a::Maybe{T}, b::Maybe{T})::Bool where {T}
    a._v == b._v
end

"""
    from_maybe_with(default::B, func, value::Maybe{T})::B where {T}

If the `Maybe` value is `nothing`, the function returns the `default` value. 
Otherwise, it applies the function to the value inside the `Maybe` and returns the result.
"""
function from_maybe_with(default::A, func, value::Maybe{T})::A where {T, A}
    is_nothing(value) ? default : func(value._v)
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
    is_nothing(value) ? error("unsafe_from_just recieved nothing value") : value._v
end

"""
    if_just_then(default, func, value::Maybe{T})

if the value is `nothing`, call the defualt function and return the result.
Otherwise, extract the value from the `Maybe`, pass it to `func`, and return the result.
"""

function if_just_then(default, func, value::Maybe{T}) where {T}
    is_nothing(value) ? default() : func(value._v)
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
    lift(func::(A -> B -> ... -> Z))::(Maybe{A} -> Maybe{B} -> ... -> Maybe{Z})

lift a function to be applicable to `Maybe` values. If any of the input values are `nothing`, `nothing` is returned.

# examples
```julia-repl
julia> maybe_add = lift(+)
(::var"#lifted#3"{typeof(+)}) (generic function with 1 method)
julia> a = maybe_add(Maybe(10), Maybe(4), Maybe(100))
Maybe{Int64}(114)
julia> unsafe_from_maybe(a)
114


julia> b = maybe_add(Maybe(10), Maybe(4), Maybe(100), EmptyMaybe)
Maybe{Nothing}(nothing)
```

"""
function lift(func)
    function lifted(vs...)::Maybe{T} where T
        any((is_nothing(v) for v in vs)) ? nothing : Maybe(func((v._v for v in vs)...))
    end
end

nothing
