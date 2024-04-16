using Pkg
Pkg.activate("..")


# for when I understand how to use Stefanos' package
# to generate the electric fields
# using ElectricFields
# using Unitful
# using UnitfulAtomic
# include("constant_envelope.jl")
# include("fields.jl")

using QuadGK
using NLsolve

include("utils.jl")
include("beams.jl")

include("equations.jl")

