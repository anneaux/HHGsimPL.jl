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

### hessian root 
function hessian_determinant(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64)
    return d2S_dtr2(b, Ip, ti, tr) * d2S_dti2(b, Ip, ti, tr)
end  

function hessian_determinant(b::Beam, Ip::Float64, s::Saddle)
    return d2S_dtr2(b, Ip, s.ti, s.tr) * d2S_dti2(b, Ip, s.ti, s.tr)
end


function dipole(b::Beam, Ip::Float64, s::Saddle) ### new 
 
  traveltime = s.tr - s.ti
 
  fontaine_SR_m0 = 1 / (kappa(Ip) * sqrt(2) * π)

  prefactor = (im * 2*π/sqrt(hessian_determinant(b, Ip, s))) # corresponds to HessianRoot in RBSFA
  prefactor *= (2*π/(im*traveltime))^(3/2) # spreading factor
  prefactor *= fontaine_SR_m0
  # transpose(conj.(dip_i)) * E(s.ti) # fontaine
  # prefactor *= sum(dip_i .* E(s.ti))

  amp = prefactor * dipole_SR_conj(s.p .+ A(b)(s.tr), Ip)

  phase = S(b, Ip, s)
  return amp .* exp(-im*phase) 
end




function harmonic_intensity(b::Beam, dipX::Complex{Float64}, dipY::Complex{Float64}, q::Number;add_cc::Bool=false)
  if add_cc
    dipX += conj(dipX)
    dipY += conj(dipY)
  end
  return (abs(dipX)^2 + abs(dipY)^2) * (q*b.omega1)^4/(2*pi*c^3)
end


function harmonic_intensity(b::Beam, Ip::Float64, s::Saddle ; add_cc::Bool=false)
  dipX, dipY = dipole(b, Ip, s)

  return harmonic_intensity(b, dipX, dipY, s.q, add_cc = add_cc)
end

