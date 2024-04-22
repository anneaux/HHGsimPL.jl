### kudos Stefanos 
using Parameters
using IntervalSets

import ElectricFields: AbstractEnvelope, envelope_types, test_field_parameters,
    continuity, span, duration, time_bandwidth_product

struct ConstantEnvelope{T} <: AbstractEnvelope
    flat::T
    period::T

    function ConstantEnvelope(flat::T, period::T) where T
        flat >= 0 || error("Negative flat region not supported")
        flat > 0 || error("Pulse length must be non-zero")

        new{T}(flat, period)
    end
end
envelope_types[:constant] = ConstantEnvelope

function (env::ConstantEnvelope{T})(t) where T
    t = real(t)/env.period
    if t < 0
        zero(T)
    elseif t <= env.flat
        one(T)
    else
        zero(T)
    end
end

show(io::IO, env::ConstantEnvelope) =
    printfmt(io, "|{1:d}| cycles Constant envelope",
             env.flat)

function ConstantEnvelope(field_params::Dict{Symbol,Any}, args...)
    test_field_parameters(field_params, [:T]) # Period time required to relate ramps/flat to cycles
    test_field_parameters(field_params, [:flat])

    @unpack flat, T = field_params
    ConstantEnvelope(flat, austrip(T))
end

continuity(::ConstantEnvelope) = 0
span(env::ConstantEnvelope{T}) where T =
    zero(T)..(env.flat*env.period)

function duration(env::ConstantEnvelope)
    s = span(env)
    s.right-s.left
end

# This is formally correct, but not very useful
time_bandwidth_product(::ConstantEnvelope) = Inf;


