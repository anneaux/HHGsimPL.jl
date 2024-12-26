


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





