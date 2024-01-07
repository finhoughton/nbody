include("utils/maybe.jl")
include("utils/partial.jl")

# ----- constants -----

const UPS::Int64 = 50

const Δt::Float64 = 1/UPS

const M_EARTH::Float64 = 6e22

const DIST::Float64 = 1e9

# const G::Float64 = 6.674e-11
const G::Float64 = 10^-3

const EPS_SOFTENING::Float64 = 1e8
# stop forcess becoming too big when objects are very close

const EDGE::Float64 = 10 * DIST

const X_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)
const Y_LIMITS::Tuple{Float64, Float64} = (-EDGE, EDGE)