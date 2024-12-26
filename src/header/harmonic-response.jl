# ### everything dipole and harmonic field


# ########

# function harmonic_field(dipole::Vector{ComplexF64}, ϕ::Float64)
#   return real.(dipole * exp(im*ϕ) )#+ conj.([dipX;dipY]) * exp(im*q*omega*t))
# end


# function harmonic_field(ellipse::Ellipse, ϕ::Float64)
#   return  real.((ellipse.A + im*ellipse.B) * exp(im * (ϕ + ellipse.theta)))
# end



# using QuadGK


kappa(Ip::Float64) = sqrt(2*Ip)
function dipole_SR_conj(k::Vector{ComplexF64}, Ip::Float64)
  ka = kappa(Ip)
  return (im *sqrt(2))/  (pi* ka) * k / (scalarproduct2(k) + ka^2 )^2
end

# ### hessian root ### F*** THIS!!! DON'T EVER DARE TO USE THIS AGAIN! IT HAS A BRANCH CUT! see Emilio's RBSFA weird handling of it
# function hessian_determinant(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64)
#   ### ffs I don't know why so far I never had this second term here
#     return d2S_dtr2(b, Ip, ti, tr) * d2S_dti2(b, Ip, ti, tr) - d2S_dtitr(b, ti, tr) * d2S_dtitr(b, ti, tr)
# end  

# function hessian_determinant(b::Beam, Ip::Float64, s::Saddle)
#     #d2S_dtr2(b, Ip, s.ti, s.tr) * d2S_dti2(b, Ip, s.ti, s.tr)
#     return hessian_determinant(b, Ip, s.ti, s.tr) 
# end
#



function hessian_root(h::AbstractArray) 
    ### I should definitely work this out properly and also make sure this actually gets rid of the branch cuts...

#     h = f_hessian(ti, tr)
    xd2S_dti2 = -conj(complex(h[1:2,1]...))
    xd2S_dtr2 = -conj(complex(h[3:4,3]...))
    xd2S_dtitr = -conj(complex(h[3:4,1]...))
    
    # ### RBSFA hessian_root
    sqrt1 = sqrt(2π/ (im*xd2S_dti2 ))
    sqrt2 = sqrt(2π * xd2S_dti2 / (im*(xd2S_dtr2 * xd2S_dti2 - xd2S_dtitr * xd2S_dtitr)) )
    return sqrt1 * sqrt2
    
    ### hessian_determinant < = works better for now. But maybe needs to be changed
    # hdet = xd2S_dtr2 * xd2S_dti2 - xd2S_dtitr * xd2S_dtitr
    return im * 2*π/sqrt(hdet)

end



function saddles_gaussian_contribution(f::Function,
#     f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64;
    prefactor::Function = (ti,tr) -> ones(2)
        )  
        
    ### prefactor for the saddle-point method
    prefactor_spm = hessian_root(f_hessian(ti,tr))

    return prefactor(ti,tr) .* prefactor_spm .* exp(f(ti,tr))
        
end 

















#### BELOW WANTS TO BE JUSTIFIED WITHIN THE NEW VERSION OF THE CODE FIRST

# function hessian_root(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64)
# #     (im * 2*π/sqrt(hessian_determinant(b, Ip, s))
#     sqrt1 = sqrt(2π/ (im*d2S_dti2(b, Ip, ti, tr) ))
#     sqrt2 = sqrt(2π * d2S_dti2(b, Ip, ti, tr) / (im*(d2S_dtr2(b, Ip, ti, tr) * d2S_dti2(b, Ip, ti, tr) - d2S_dtitr(b, ti, tr) * d2S_dtitr(b, ti, tr))) )
    
#     return sqrt1 * sqrt2       
        
# end

# function hessian_root(b::Beam, Ip::Float64, s::Saddle)
#     return hessian_root(b, Ip, s.ti, s.tr)
# end



# function dipole(b::Beam, Ip::Float64, s::Saddle) ### new 
 
#   traveltime = s.tr - s.ti
 
#   fontaine_SR_m0 = 1 / (kappa(Ip) * sqrt(2) * π)

#   prefactor = hessian_root(b, Ip, s) # corresponds to HessianRoot in RBSFA
#   prefactor *= (2*π/(im*traveltime))^(3/2) # spreading factor
#   prefactor *= fontaine_SR_m0

#   # transpose(conj.(dip_i)) * E(s.ti) # fontaine
#   # prefactor *= sum(dip_i .* E(s.ti))

#   amp = prefactor * dipole_SR_conj(s.p .+ A(b)(s.tr), Ip)

#   phase = S(b, Ip, s)
#   return amp .* exp(-im*phase) 
# end


# function harmonic_intensity(b::Beam, dipX::Complex{Float64}, dipY::Complex{Float64}, q::Number;add_cc::Bool=false)
#   if add_cc
#     dipX += conj(dipX)
#     dipY += conj(dipY)
#   end
#   return (abs(dipX)^2 + abs(dipY)^2) * (q*b.omega1)^4/(2*pi*c^3)
# end


# function harmonic_intensity(b::Beam, Ip::Float64, s::Saddle ; add_cc::Bool=false)
#   dipX, dipY = dipole(b, Ip, s)

#   return harmonic_intensity(b, dipX, dipY, s.q, add_cc = add_cc)
# end

nothing