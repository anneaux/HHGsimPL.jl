using Pkg
Pkg.activate("..")

using QuadGK
using Sobol
using NLsolve
using Plots

using LinearAlgebra
using Integrals
using FiniteDiff
using StaticArrays
using Contour




include("utils.jl")
include("beams.jl")

include("saddles.jl")
include("equations.jl")

include("stokes-and-cutoffs.jl")
include("necklace.jl")
include("saddles-contributing.jl")
