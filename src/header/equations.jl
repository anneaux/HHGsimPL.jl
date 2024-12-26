## stationary momentum
function p_stationary(b::Beam, 
    ti::ComplexF64, tr::ComplexF64)

    -1/(tr-ti) * IA(b)(ti,tr)
end

### a comment on the square here to resolve my insecurities once and for all:
# what I want is a 'fake' scalar product, that takes the total of the element-wise product of the two vectors. I think this is because originally we weren't thinking of it to be in the complex plane, hence we don't use a proper dot product (which would be sum(conj(zvec) .* zvec)) ).

## Volkov action
function S_v(b::Beam, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64, 
  p::Vector{ComplexF64} = p_stationary(b, ti, tr)
  )

  0.5 * quadgk(t -> scalarproduct2( p .+ A(b)(t) ), ti, tr)[1] + Ip * (tr - ti)
end

## exponent = action
# function S(b::Beam, Ip::Float64, 
#   ti::ComplexF64, tr::ComplexF64, 
#   q::Number,
#   p::Vector{ComplexF64} = p_stationary(b, ti, tr)
#   )
  
#   S_v(b, Ip, ti, tr , p) - q * b.omega1 * tr

# end 

function S(b::Beam, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64, 
  q::Number,
  p::Vector{ComplexF64} = p_stationary(b, ti, tr)
  )
    
    try
        S_v_analytic(b, Ip, ti, tr) - q * b.omega1 * tr
    catch e
        S_v(b,Ip,ti,tr,p) - q * b.omega1 * tr
    end
end


S_v(b::Beam, Ip::Float64, s::Saddle) = S_v(b, Ip, s.ti, s.tr, s.p) 
S(b::Beam, Ip::Float64, s::Saddle) = S(b, Ip, s.ti, s.tr, s.q, s.p)


### derivatives of the Volkov action
function dSv_dtr(b::Beam, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64) 
  # inv_tau = 1/(tr - ti)
  # IA = integral_over_A(ti,tr)
  kr = p_stationary(b, ti, tr) .+ A(b)(tr)

  return 0.5 * scalarproduct2(kr) + Ip
end 

function dSv_dti(b::Beam, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64)
  # inv_tau = 1/(tr - ti)
  # IA = integral_over_A(ti,tr)
  ki = p_stationary(b, ti, tr) .+ A(b)(ti)

  return  -1 * (0.5 * scalarproduct2(ki) + Ip)
end

### derivatives of the action
function dS_dtr(b::Beam, Ip::Float64,
  q::Number, 
  ti::ComplexF64,tr::ComplexF64) 

  return dSv_dtr(b, Ip, ti, tr) - q* b.omega1
end 
 
function dS_dti(b::Beam, Ip::Float64,
  ti::ComplexF64,tr::ComplexF64) 

  return  dSv_dti(b, Ip, ti, tr)
end

# just for consistency of definitions
function dS_dti(b::Beam, Ip::Float64,
  q::Number, 
  ti::ComplexF64, tr::ComplexF64) 
  return  dSv_dti(b, Ip, ti, tr)
end



### second drv of the Volkov action
function d2Sv_dtr2(b::Beam,
  ti::ComplexF64,tr::ComplexF64)

  inv_tau = 1/(tr - ti)
  kr = p_stationary(b, ti, tr) .+ A(b)(tr)

  return scalarproduct( kr , - E(b)(tr) .- inv_tau .*  kr)
end


function d2Sv_dti2(b::Beam, 
  ti::ComplexF64,tr::ComplexF64) 

  ki = p_stationary(b, ti, tr) .+ A(b)(ti)
  inv_tau = 1/(tr - ti)

  return scalarproduct( ki , E(b)(ti) .- inv_tau .* ki)
end 


### mixed second derivative
function d2Sv_dtitr(b::Beam, 
  ti::ComplexF64, tr::ComplexF64) 
  #pfun = -1/(tr-ti) * integral_over_A(ti,tr)
  inv_tau = 1/(tr - ti)
  integral = IA(b)(ti,tr)

  return inv_tau * 
    ( inv_tau * scalarproduct(- (A(b)(ti) .+ A(b)(tr)), integral ) + 
    + scalarproduct( A(b)(ti), A(b)(tr) ) + 
    + inv_tau *inv_tau * scalarproduct2(integral)
    )
end

d2S_dtitr(b::Beam, ti::ComplexF64, tr::ComplexF64) = d2Sv_dtitr(b, ti, tr)
d2S_dtitr(b::Beam, s::Saddle) = d2Sv_dtitr(b, s.ti, s.tr)



### just for consistency in the code
  d2Sv_dti2(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2Sv_dti2(b::Beam, ti::ComplexF64, tr::ComplexF64)

  d2Sv_dtr2(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2Sv_dtr2(b::Beam, ti::ComplexF64,tr::ComplexF64) 

  d2Sv_dtitr(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2Sv_dtitr(b::Beam, ti::ComplexF64,tr::ComplexF64)
  d2S_dtitr(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2S_dtitr(b::Beam, ti::ComplexF64, tr::ComplexF64)


### second derivatives of action S
  d2S_dti2(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2Sv_dti2(b::Beam, ti::ComplexF64, tr::ComplexF64)
  d2S_dti2(b::Beam, ti::ComplexF64, tr::ComplexF64) = d2Sv_dti2(b::Beam, ti::ComplexF64, tr::ComplexF64)

  d2S_dtr2(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = d2Sv_dtr2(b::Beam, ti::ComplexF64, tr::ComplexF64)
  d2S_dtr2(b::Beam, Ip::Float64, q::Number, ti::ComplexF64, tr::ComplexF64) = d2S_dtr2(b, Ip, ti, tr)


### third drv of the action
# TODO 

#### equations for the necklace

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























## Saddle point equations (SPEQs)

#### I DON'T THINK I NEED THEM ANYMORE
# this is actually just S_V_drv(p,ti)
# function speq1(b::Beam, Ip::Float64,
#   tir::Float64, tii::Float64, trr::Float64, tri::Float64
#   )
#   ti = tir + im * tii
#   tr = trr + im * tri
  
#   0.5 * scalarproduct2( p_stationary(b, ti, tr) + A(b)(ti) ) .+ Ip

# end 

# speq1(b::Beam, Ip::Float64, ti::ComplexF64, tr::ComplexF64) = speq1(b, Ip, real(ti), imag(ti),real(tr), imag(tr))

# # this is actually just S_drv = S_V_drv(p,tr) - q*omega
# function speq2(b::Beam, Ip::Float64,
#   q::Number, 
#   tir::Float64, tii::Float64, trr::Float64, tri::Float64
#   )
#   ti = tir + im * tii
#   tr = trr + im * tri  

#   0.5 * scalarproduct2( p_stationary(b, ti, tr) + A(b)(tr) ) .+ Ip .- q * b.omega1
# end 

# speq2(b::Beam, Ip::Float64, q::Number, ti::ComplexF64, tr::ComplexF64) = speq2(b, Ip, q, real(ti),imag(ti),real(tr),imag(tr))


  nothing