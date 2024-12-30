struct BeamMono <: Beam
	E0::Float64
	omega1::Float64
	epsilon::Float64
	phi::Float64

	function BeamMono(;Intensity::Float64,
	 lambda::Int64,
	 epsilon::Float64=0.,
	 phi::Float64=0.)
	
		E0 = sqrt(Intensity/IAU)
		omega1 = get_omega(lambda)
		
		new(E0,omega1,epsilon,phi)
	end


	function BeamMono(E0::Float64,
		omega1::Float64,
		epsilon::Float64,
		phi::Float64)

		new(E0,omega1,epsilon,phi)
	end

end

TCycle(b::BeamMono) = 2*pi/b.omega1

function electric_field(b::BeamMono)
	E(t) = (b.E0/sqrt(1+b.epsilon^2) .* [sin(b.omega1*t+ b.phi) ; -b.epsilon*cos(b.omega1*t+ b.phi)])
	return t -> E(t)
end
E(b::BeamMono) = electric_field(b)


function vector_potential(b::BeamMono)
	A0 = b.E0/b.omega1/sqrt(1+b.epsilon^2) 
	A(t) = (A0 .* [cos(b.omega1*t + b.phi) ; b.epsilon * sin(b.omega1*t + b.phi)])
	
	return t -> A(t)
end
A(b::BeamMono) = vector_potential(b)



function integrated_vector_potential(b::BeamMono)
	A0 = b.E0/b.omega1/sqrt(1+b.epsilon^2)

	integral_over_A(ti,tr) = (A0 ./ b.omega1 .* [sin(b.omega1*tr+ b.phi) - sin(b.omega1*ti+ b.phi); b.epsilon*(-cos(b.omega1*tr+ b.phi) + cos(b.omega1 * ti+ b.phi))])
	
	return (ti, tr) -> integral_over_A(ti,tr)
end
IA(b::BeamMono) = integrated_vector_potential(b)


function integrated_squared_vector_potential(b::BeamMono)
    A0 = b.E0/b.omega1/sqrt(1+b.epsilon^2)
    A0sq = A0^2
    omega = b.omega1
    phi = b.phi
    
    (ti,tr) ->  A0sq/(4*omega) * (2*(1+epsilon^2) * omega * (tr - ti) + (-1+epsilon^2) * sin(2*(omega*ti + phi)) - (-1+epsilon^2)*sin(2*(omega*tr + phi)))

end
IAsq(b::BeamMono) = integrated_squared_vector_potential(b)


function S_v_analytic(b::BeamMono, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64, 
  ps::Vector{ComplexF64} = p_stationary(b, ti, tr)
  )
    pssq = scalarproduct2(ps)
#     integral = integrate( sum((p+A(t)) .* (p+A(t)))  , t, ti, tr)
    
    integral = pssq*(tr-ti) + 2 * sum(ps .* IA(b)(ti,tr)) + IAsq(b)(ti,tr)
   
    0.5 * integral  + Ip * (tr - ti)
end

S_v_analytic(b::BeamMono, Ip::Float64, s::Saddle) = S_v_analytic(b, Ip, s.ti, s.tr, s.p);

# function integral_over_E(b::BeamMono,t::Float64)
# 	res= b.E0/sqrt(1+b.epsilon^2) .*  [-((cos(b.phi)*cos(b.omega1*t))/b.omega1)+(sin(b.phi)*sin(b.omega1*t))/b.omega1 ; -b.epsilon*((cos(b.omega1*t)*sin(b.phi))/b.omega1+(cos(b.phi)*sin(b.omega1*t))/b.omega1)]
# 	return res
# end

# function second_integral_over_E(b::BeamMono,t::Float64)
# 	res = b.E0/sqrt(1+b.epsilon^2) .* [-(((cos(b.omega1*t)*sin(b.phi))/b.omega1+(cos(b.phi)*sin(b.omega1*t))/b.omega1)/b.omega1) ; -((b.epsilon*(-((cos(b.phi)*cos(b.omega1*t))/b.omega1)+(sin(b.phi)*sin(b.omega1*t))/b.omega1))/b.omega1)]

# 	return res
# end

function get_Up(b::BeamMono)
	A0 = b.E0/b.omega1
	Up = A0^2/4 
	return Up
end


