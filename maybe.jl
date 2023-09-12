
struct Maybe{T}
    _v::Union{T, Nothing}
end

const EmptyMaybe = Maybe(nothing)

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
    conditional_to_maybe(default::T, func, value::T)::Maybe{T} where {T}

if the function applied to the values returns `True`,
a `Maybe` value containing the value is returned, otherwise a `Maybe` value containing the defualt is returned
"""
function conditional_to_maybe(default::T, func, value::T)::Maybe{T} where {T}
    func(value) ? Maybe(value) : Maybe(default)
end

"""
    unsafe_from_maybe(value::Maybe{T})::T where {T}

Attempts to extract the value from a `Maybe`.
if it is `nothing`, an error is raised, otherwise the value is returned.
"""
function unsafe_from_maybe(value::Maybe{T})::T where {T}
    is_nothing(o) ? error("unsafe_from_maybe recieved nothing value") : value._v
end

"""
    if_maybe_then(default, func, value::Maybe{T})

if the value is `nothing`, call the defualt function and return the result.
Otherwise, extract the value from the `Maybe`, pass it to `func`, and return the result.
"""

function if_maybe_then(default, func, value::Maybe{T}) where {T}
    is_nothing(value) ? default() : func(value._v)
end

"""
    is_something(value::Maybe{T})::Bool where {T}

if a `Maybe` value contains a value
"""
function is_something(value::Maybe{T})::Bool where {T}
    value._v ≢ nothing
end

"""
    is_nothing(value::Maybe{T})::Bool where {T}

if a `Maybe` value is `nothing`
"""
function is_nothing(value::Maybe{T})::Bool where {T}
    value._v ≡ nothing
end

"""
    maybe_map(func, value::Maybe{T})::Maybe{S} where {T, S}

Apply a function under a `Maybe` value, if the value is nothing, it is unchanged.
"""
function maybe_map(func, value::Maybe{T})::Maybe{T} where {T}
    (Maybe ∘ from_maybe_with)(nothing, func, value)
end
