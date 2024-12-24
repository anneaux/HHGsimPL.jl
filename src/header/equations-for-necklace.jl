# I guess I could at some point use this for the usual calculation of SV as well
function S_v_for_diff(b::Beam, Ip::Float64, 
  ti::Complex, tr::Complex, 
  p::Vector = p_stationary(b, ti, tr)
  )
    domain = (ti,tr)
    prob = IntegralProblem((t,x) -> scalarproduct2( p .+ A(b)(t) ), domain)
    integral = solve(prob, QuadGKJL())[1]
  return 0.5 * integral + Ip * (tr - ti)
end

function S_for_diff(b::Beam, Ip::Float64, 
  ti::Complex, tr::Complex, 
  q::Number,
  p::Vector = p_stationary(b, ti, tr)
  )
  
  S_v_for_diff(b, Ip, ti, tr , p) - q* b.omega1 * tr
end


function my_hessian(b::Beam, Ip::Float64,
        q::Number,
        ti::ComplexF64, tr::ComplexF64)

    action(tvec) = real(-S_for_diff(b, Ip, tvec[1]+im*tvec[2], tvec[3]+im*tvec[4], q)) # doesn't matter if I take real or imag there
    return FiniteDiff.finite_difference_hessian(action, [reim(ti)..., reim(tr)...])
end



### grads and hessians

	# function grad(b::Beam, Ip::Float64,
	#     q::Number,
	#     ti::ComplexF64, tr::ComplexF64)

	#     g = [dS_dti(b,Ip,q,ti,tr); dS_dtr(b,Ip,q,ti,tr)]
	#     g = conj.(complex.(-1im .* g))
	#     return g
	# end

########integrate quadrilaterals
	
# ti,tr = map([x[i], x[j]], p1, p2, p3, p4)
#             action = f_vec([ti,tr])
            
#             ### eq. 5 from Emilio's thc paper
#             ps = p_stationary(b, ti, tr)
#             dip_r = dipole_SR_conj(ps .+ A(b)(tr), Ip)
#             fontaine_SR_m0 = 1 / (kappa(Ip) * sqrt(2) * π)
#             traveltime = tr - ti
            
#             prefactor = (2*π/(im*traveltime))^(3/2) # spreading factor
#             prefactor *= fontaine_SR_m0

#             sum = sum + jac * prefactor * dip_r * exp(-im * action) * w[i] * w[j]





