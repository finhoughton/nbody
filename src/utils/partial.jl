

"""
partial function application in julia, same as python's functools.partial
"""
function partial(f, a...)
    ((b...) -> f(a..., b...))
end
