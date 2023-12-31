

"""
partial function application in julia, same as python's functools.partial
"""
function partial(f, a...)
    function wrapped(b...)
        return f(a..., b...)
    end
    return wrapped
end
