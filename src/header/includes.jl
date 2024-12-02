using Pkg
Pkg.activate("..")

using QuadGK
using Sobol
using NLsolve
# using Plots

using LinearAlgebra
using Integrals
using FiniteDiff
using StaticArrays
using Contour
using Peaks



include("utils.jl")
include("beams.jl")
include("Beams/BeamMono.jl")
include("Beams/BeamTC.jl")
include("Beams/BeamOTC.jl")


include("saddles.jl")
include("equations.jl")

include("stokes-and-cutoffs.jl")
include("necklace.jl")
include("saddles-contributing.jl")
include("harmonic-response.jl")
