using Pkg
Pkg.activate("..")



using Sobol
using NLsolve
using QuadGK
using Integrals
using FiniteDiff
using StaticArrays
using Contour
using Peaks
using LinearAlgebra
using GeometryBasics


include("constants-units.jl")
include("utils.jl")

include("beams.jl")
include("saddles.jl")

include("Beams/BeamMono.jl")
include("Beams/BeamTC.jl")
include("Beams/BeamOTC.jl")


include("equations.jl")
include("harmonic-response.jl")
include("stokes-and-cutoffs.jl")
include("ionisationbursts.jl")
