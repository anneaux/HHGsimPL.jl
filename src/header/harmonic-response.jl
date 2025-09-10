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

function dipole_SR_conj(k::AbstractVector{ComplexF64}, Ip::Float64)
  ka = kappa(Ip)
  return (im *sqrt(2))/  (pi* ka) * k / (scalarproduct2(k) + ka^2 )^2
end




function hessian_root(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64)
#     (im * 2*π/sqrt(hessian_determinant(b, Ip, s))
    sqrt1 = sqrt(2π/ (im*d2S_dti2(b, Ip, ti, tr) ))
    sqrt2 = sqrt(2π * d2S_dti2(b, Ip, ti, tr) / (im*(d2S_dtr2(b, Ip, ti, tr) * d2S_dti2(b, Ip, ti, tr) - d2S_dtitr(b, ti, tr) * d2S_dtitr(b, ti, tr))) )
    
    return sqrt1 * sqrt2            
end

function hessian_root(b::Beam, Ip::Float64, s::Saddle)
    return hessian_root(b, Ip, s.ti, s.tr)
end


function dipole(b::Beam, Ip::Float64, s::Saddle) ### new 
 
  traveltime = s.tr - s.ti
 
  fontaine_SR_m0 = 1 / (kappa(Ip) * sqrt(2) * π)

  prefactor = hessian_root(b, Ip, s) # corresponds to HessianRoot in RBSFA
  prefactor *= (2*π/(im*traveltime))^(3/2) # spreading factor
  prefactor *= fontaine_SR_m0

  # transpose(conj.(dip_i)) * E(s.ti) # fontaine
  # prefactor *= sum(dip_i .* E(s.ti))

  amp = prefactor * dipole_SR_conj(s.p .+ A(b)(s.tr), Ip)

  phase = S(b, Ip, s)
  return amp .* exp(-im*phase) 
end



#####


function harmonic_intensity(b::Beam, dip::AbstractVector{ComplexF64}, q::Number; add_cc::Bool=false)
  if add_cc
    dip[1] += conj(dip[1])
    dip[2] += conj(dip[2])
  end
  return sum((abs.(dip)).^2) * (q*fundamental_frequency(b))^4/(2*pi*c^3)
end

nothing