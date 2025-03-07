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


function harmonic_intensity(b::Beam, dip::AbstractVector{ComplexF64}, q::Number; add_cc::Bool=false)
  if add_cc
    dip[1] += conj(dip[1])
    dip[2] += conj(dip[2])
  end
  return sum((abs.(dip)).^2) * (q*fundamental_frequency(b))^4/(2*pi*c^3)
end

nothing