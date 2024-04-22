include("../src/header/includes.jl")
using Plots
using LaTeXStrings


# @field(IR) do
#     I₀ = 1e14u"W/cm^2"
#     λ = 800.0u"nm"
#     τ = 6.2u"fs"
#     σmax = 6.0
# end

field = get_field_XZ(intratio=2.,ϕ=0.)

au2fs = ustrip(auconvert(u"fs", 1.0))
t = timeaxis(field)
tplot = au2fs*t

Fv = field_amplitude(field, t)
Av = vector_potential(field, t)

Fp = plot(tplot, Fv, ylabel=L"$F(t)$ [au]")
Ap = plot(tplot, Av, ylabel=L"$A(t)$ [au]")
plot(Fp, Ap, layout=@layout([a;b]), xlabel=L"$t$ [fs]")

