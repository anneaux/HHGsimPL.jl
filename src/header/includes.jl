using Pkg
Pkg.activate("..")

using QuadGK
using NLsolve


include("utils.jl")
include("beams.jl")

include("equations.jl")

include("saddles.jl")
