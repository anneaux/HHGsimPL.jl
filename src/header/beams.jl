abstract type Beam end

TCycle(b::Beam) = 2*pi/fundamental_frequency(b)

# TC(lambda::Int64) = 2*pi/get_omega(lambda)
TCycle(;lambda::Real) = lambda / (LAU * c)
TCycleNU(;lambda::Real) = lambda * 1e-9/cNU
get_omega(lambda::Real) = 2 * pi * LAU * c / lambda


function field_amplitude(vec::Vector{T}) where T <: Real
	return sqrt(sum(vec.^2))
end


