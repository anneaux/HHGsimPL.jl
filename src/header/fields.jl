### here I can collect some common field definitions.
# e.g. switchover, OTC, BEOTC, BICIRC, MONO

# function get_field_XZ(;intratio::Number=0,phi::Number=0)
#     @field(A) do
#         I₀ = 1e14u"W/cm^2"
#         λ = 800.0u"nm"
#         ramp = 0.0
#         flat = 1.0
#         env = :trapezoidal
#     end
#     @field(B) do
#         I₀ = intratio*1e14u"W/cm^2"
#         λ = 400.0u"nm"
#         ramp = 0.0
#         flat = 2.0
#         env = :trapezoidal
#         ϕ = phi
#         rotation = π/2, [0,0,1]
#     end

# 	A + B
# end

### this is a quick hack until Stefanos tells me how to do it properly
import ElectricFields: AbstractField

function integrated_vector_potential(f::AbstractField,t1::Number, t2::Number)
	quadgk(t -> vector_potential(f,t) ), t1, t2)[1]
end 



function get_field_XZ(;intratio::Number=1, ϕ::Number=0)
    @field(A) do
        I₀ = 1e14u"W/cm^2"
        λ = 800.0u"nm"
        flat = 1.0
        env = :constant
    end
    @field(B) do
        I₀ = intratio*1e14u"W/cm^2"
        λ = 400.0u"nm"
        flat = 2.0
        env = :constant
        ф = ϕ
        rotation = π/2, [1,0,0]
    end
   
   A + B
end


# function get_field_TC(;theta::Number=0,phi::Number=0)
#     s,c = sincos(theta)
#     @field(A) do
#         I₀ = c*1e14u"W/cm^2"
#         λ = 800.0u"nm"
#         ramp = 0.0
#         flat = 1.0
#         env = :trapezoidal
#     end
#     @field(B) do
#         I₀ = s*1e14u"W/cm^2"
#         λ = 400.0u"nm"
#         ramp = 0.0
#         flat = 2.0
#         env = :trapezoidal
#         ϕ = phi
#     end

# 	A + B
# end

# @field(A) do
#     I₀ = c*1e14u"W/cm^2"
#     λ = 800.0u"nm"
#     τ = 6.2u"fs"
#     σmax = 6.0
# end
# @field(B) do
#     I₀ = s*1e14u"W/cm^2"
#     λ = 400.0u"nm"
#     τ = 6.2u"fs"
#     σmax = 6.0
# end
# @field(A) do
#     I₀ = c*1e14u"W/cm^2"
#     λ = 800.0u"nm"
#     τ = 6.2u"fs"
#     σoff = 4.0
#     σmax = 6.0
#     env = :trunc_gauss
# end
# @field(B) do
#     I₀ = s*1e14u"W/cm^2"
#     λ = 400.0u"nm"
#     τ = 6.2u"fs"
#     σoff = 4.0
#     σmax = 6.0
#     env = :trunc_gauss
# end

# A + phase_shift(B, phi)



### integral over the vector potential