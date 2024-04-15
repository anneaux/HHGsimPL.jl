# ### everything related to saddle points

# ### find saddles
# ### different methods: classical guess, sobol, nroots



# abstract type SP end

# ### Saddle ###########################
# struct Saddle <: SP
# 	q::T where T <: Number
# 	ti::Complex{Float64}
# 	tr::Complex{Float64}
# 	p::Vector{ComplexF64} # do I want this here???

# end

# Base.show(io::IO, s::Saddle) = print(io,
#     "ti: $(round(s.ti,sigdigits=5)), tr: $(round(s.tr,sigdigits=5)), q$(round(s.q,sigdigits=5))")

# import Base./
# function /(s::Saddle,TC::Real)
#     return Saddle(s.q, s.ti/TC, s.tr/TC, s.p)
# end
