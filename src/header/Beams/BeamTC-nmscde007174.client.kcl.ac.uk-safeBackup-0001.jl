#### struct and constructors ###########
struct BeamTC <: Beam
	E01::Float64
	E02::Float64
	omega1::Float64
	omega2::Float64
	epsilon1::Float64
	epsilon2::Float64
	phi::Float64

	function BeamTC(;Intensity1::Real, Intensity2::Real,
		# omega
		lambda::Real, 
		r::Int64, s::Int64,
		epsilon1::Real, epsilon2::Real,
		phi::Real)
		### all below is just overflow
		# get_omega(lambda)

		E01 = sqrt(Intensity1/IAU)
		E02 = sqrt(Intensity2/IAU)
		omega = get_omega(lambda)
		omega1 = r*omega
		omega2 = s*omega
		epsilon1 = epsilon1
		epsilon2 = epsilon2
		
		new(E01,E02,omega1,omega2,epsilon1,epsilon2,phi)
	end

	function BeamTC(E01::Float64,
		E02::Float64,
		omega1::Float64,
		omega2::Float64,
		epsilon1::Float64,
		epsilon2::Float64,
		phi::Float64)
		
		new(E01,E02,omega1,omega2,epsilon1,epsilon2,phi)
	end
end

# outer constructor for equal-ellipticities
function BeamBETC(;Intensity1::Real, Intensity2::Real,
	lambda::Real,
	r::Int64,s::Int64,
	epsilon::Real,
	phi::Real)

	BeamTC(Intensity1=Intensity1,Intensity2=Intensity2,
		lambda=lambda,
		r=r,s=s,
		epsilon1=epsilon,epsilon2=epsilon,
		phi=phi)
end

TCycle(b::BeamTC) = 2*pi/b.omega1


#### field equations ###########


function electric_field(b::BeamTC)
	## Beam Characterization from milo2020biell
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	rotmat = 1 #[0 -1; -1 0]

	E1(t) = rotmat * (E01/sqrt(1+epsilon1^2) .* [sin(omega1*t) ; -epsilon1*cos(omega1*t)] )
	E2(t) = rotmat * (E02/sqrt(1+epsilon2^2) .* [sin(omega2*t + phi); -epsilon2*cos(omega2*t + phi)] )

	E(t) = E1(t) + E2(t)

	return t -> E(t)
end
E(b::BeamTC) = electric_field(b)


function vector_potential(b::BeamTC)
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	rotmat = 1 #[0 -1; -1 0]

	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)
	A1(t) = rotmat * (A01 .* [cos(omega1*t) ; epsilon1 * sin(omega1*t)])
	A2(t) = rotmat * (A02 .* [cos(omega2*t + phi); epsilon2 * sin(omega2*t + phi) ])

	A(t) = A1(t) .+ A2(t)

	return t -> A(t)
end 
A(b::BeamTC) = vector_potential(b)



function integrated_vector_potential(b::BeamTC)

	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	rotmat = 1 #[0 -1; -1 0]

	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)

	integral_over_A1(ti,tr) = rotmat * (A01 ./ omega1 .* [sin(omega1*tr) - sin(omega1*ti); epsilon1*(-cos(omega1*tr) + cos(omega1 * ti))])
	integral_over_A2(ti,tr) = rotmat * (A02 ./ omega2 .* [ sin(omega2*tr + phi) - sin(omega2*ti + phi) ; epsilon2*(-cos(omega2*tr + phi) + cos(omega2*ti + phi))])

	integral_over_A(ti,tr) = integral_over_A1(ti,tr) .+ integral_over_A2(ti,tr)

	return (ti,tr) -> integral_over_A(ti,tr)
end 

IA(b::BeamTC) = integrated_vector_potential(b)


function integrated_squared_vector_potential(b::BeamTC)
    phi = b.phi
    epsilon1 = b.epsilon1
    E01 = b.E01
    omega1 = b.omega1
    E02 = b.E02
    epsilon2 = b.epsilon2
    omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

    A01 = E01/omega1/sqrt(1+epsilon1^2) 
    A02 = E02/omega2/sqrt(1+epsilon2^2)
    
    A01sq = A01^2
    A02sq = A02^2
    
    term1(ti,tr) = A02^2 * (-1 + epsilon2^2) * omega1 * (omega1 - omega2) * (omega1 + omega2) * (cos(2 * omega2 * ti) - cos(2 * omega2 * tr)) * sin(2 * phi)
    term2(ti,tr) = 8 * A01 * A02 * omega1 * omega2 * (epsilon1 * epsilon2 * omega1 + omega2) * cos(omega1 * ti) * sin(phi + omega2 * ti)
    term3(ti,tr) = -8 * A01 * A02 * omega1 * omega2 * (omega1 + epsilon1 * epsilon2 * omega2) * cos(phi) * (cos(omega2 * ti) * sin(omega1 * ti) - cos(omega2 * tr) * sin(omega1 * tr))
    term4(ti,tr) = omega2 * (-omega1^2 + omega2^2) * (2 * (A01^2 * (1 + epsilon1^2) + A02^2 * (1 + epsilon2^2)) * omega1 * (ti - tr) + A01^2 * (-1 + epsilon1^2) * (-sin(2 * omega1 * ti) + sin(2 * omega1 * tr)))
    term5(ti,tr) = 8 * A01 * A02 * omega1 * omega2 * (omega1 + epsilon1 * epsilon2 * omega2) * sin(phi) * (sin(omega1 * ti) * sin(omega2 * ti) - sin(omega1 * tr) * sin(omega2 * tr))
    term6(ti,tr) = A02^2 * (-1 + epsilon2^2) * omega1 * (omega1 - omega2) * (omega1 + omega2) * cos(2 * phi) * (sin(2 * omega2 * ti) - sin(2 * omega2 * tr))
    term7(ti,tr) = -8 * A01 * A02 * omega1 * omega2 * (epsilon1 * epsilon2 * omega1 + omega2) * cos(omega1 * tr) * sin(phi + omega2 * tr)

    (ti,tr) -> (term1(ti,tr) + term2(ti,tr) + term3(ti,tr) + term4(ti,tr) + term5(ti,tr) + term6(ti,tr) + term7(ti,tr)) / (4 * omega1 * (omega1 - omega2) * omega2 * (omega1 + omega2))

end
IAsq(b::BeamTC) = integrated_squared_vector_potential(b)


function S_v_analytic(b::BeamTC, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64, 
  ps::Vector{ComplexF64} = p_stationary(b, ti, tr)
  )
    pssq = scalarproduct2(ps)
#     integral = integrate( sum((p+A(t)) .* (p+A(t)))  , t, ti, tr)
    integral = pssq*(tr-ti) + 2 * sum(ps .* IA(b)(ti,tr)) + IAsq(b)(ti,tr)
    
    0.5 * integral  + Ip * (tr - ti)
end

######################
### indefinite integrals (required for ATI)

function integrated_vector_potential_indefinite(b::BeamTC)

	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	rotmat = 1 #[0 -1; -1 0]

	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)

	integral_over_A1(t) = rotmat * (A01 ./ omega1 .* [sin(omega1*t); epsilon1*(-cos(omega1*t))])
	integral_over_A2(t) = rotmat * (A02 ./ omega2 .* [ sin(omega2*t + phi) ; epsilon2*(-cos(omega2*t + phi) )])

	integral_over_A(t) = integral_over_A1(t) .+ integral_over_A2(t)

	return t -> integral_over_A(t)
end 

IA_indefinite(b::BeamTC) = integrated_vector_potential_indefinite(b)



function integrated_squared_vector_potential_indefinite(b::BeamTC)
    phi = b.phi
    epsilon1 = b.epsilon1
    E01 = b.E01
    omega1 = b.omega1
    E02 = b.E02
    epsilon2 = b.epsilon2
    omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

    A01 = E01/omega1/sqrt(1+epsilon1^2) 
    A02 = E02/omega2/sqrt(1+epsilon2^2)

		term1(t) = 2 * t * (A01^2 * (1 + epsilon1^2) + A02^2 * (1 + epsilon2^2)) 
		term2(t) = - (A02^2 * (-1 + epsilon2^2) * sin(2 * (phi + omega2 * t))) / omega2
		term3(t) = (A01 * (-1 + epsilon1^2) * sin(2 * omega1 * t)) / omega1
		term4(t) = (4 * A02 * (1 + epsilon1 * epsilon2) * sin(phi - omega1 * t + omega2 * t)) / (omega1 - omega2)
		term5(t) = (4 * A02 * (-1 + epsilon1 * epsilon2) * sin(phi + (omega1 + omega2) * t)) / (omega1 + omega2)

	return t -> (1 / 4) * (term1(t) + term2(t) - A01 * (term3(t) + term4(t) + term5(t)))


end
IAsq_indefinite(b::BeamTC) = integrated_squared_vector_potential_indefinite(b)




### needs to be reworked for TC beam!!!
# function electric_field_amplitude_derivative(b::BeamOTC)
# 	## Beam Characterization from milo2020biell
# 	phi = b.phi
# 	epsilon1 = b.epsilon1
# 	E01 = b.E01
# 	omega1 = b.omega1
# 	E02 = b.E02
# 	epsilon2 = b.epsilon2
# 	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

# 	term1 = 1 + epsilon1^2
# 	term2 = 1 + epsilon2^2

#    numerator(t) = (-E01 * E02 * sqrt(term1) * sqrt(term2) * (epsilon2 * omega1 + epsilon1 * omega2) * cos(t * omega1) * cos(phi + t * omega2) -
#    	0.5 * E01^2 * (-term1) * term2 * omega1 * sin(2 * t * omega1) -
#    	E02 * (E02 * term1 * (-term2) * omega2 * cos(phi + t * omega2) - 
#    	E01 * sqrt(term1) * sqrt(term2) * (epsilon1 * omega1 + epsilon2 * omega2) * sin(t * omega1))   	* sin(phi + t * omega2))

#    denominator(t) = term1 * term2 * sqrt(((E02 * epsilon2 * cos(phi + t * omega2)) / sqrt(term2) - (E01 * sin(t * omega1)) / 	sqrt(term1))^2 + 
#    	((E01 * epsilon1 * cos(t * omega1)) / sqrt(term1) - (E02 * sin(phi + t * omega2)) / sqrt(term2))^2)

#     return t -> (numerator(t) / denominator(t))
# end



function get_Up(b::BeamTC)
	A01 = b.E01/b.omega1
	A02 = b.E02/b.omega2
	# println("Up = $(A01^2/4 + A02^2/4)")
	Up = A01^2/4 + A02^2/4
	return Up
end
