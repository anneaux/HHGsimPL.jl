# Complex domain

import Base.+, Base.*, Base.in

struct ComplexDomain
    min::ComplexF64
    max::ComplexF64 #Union{ComplexF64,Nothing}
    
    ComplexDomain(rmin::Real,rmax::Real,imin::Real,imax::Real) = new(rmin+imin*im,rmax+imax*im)
    
    ComplexDomain(min::ComplexF64,max::ComplexF64) = new(min,max)

    ComplexDomain() = new(zero(ComplexF64),zero(ComplexF64))
    
end

function +(cd1::ComplexDomain,cd2::ComplexDomain)
    rmin = minimum([real(cd1.min),real(cd2.min)])
    rmax = maximum([real(cd1.max),real(cd2.max)])
    imin = minimum([imag(cd1.min),imag(cd2.min)])
    imax = maximum([imag(cd1.max),imag(cd2.max)])
    
    return ComplexDomain(rmin,rmax,imin,imax)
end

function +(cd1::ComplexDomain,z::Number) # this function shifts the whole domain by the number specified
    
    rmin = real(cd1.min) + real(z)
    rmax = real(cd1.max) + real(z)
    imin = imag(cd1.min) + imag(z)
    imax = imag(cd1.max) + imag(z)
    
    return ComplexDomain(rmin,rmax,imin,imax)
end

function *(cd1::ComplexDomain,z::Real) # this function shifts the whole domain by the number specified
    
    rmin = real(cd1.min) * z
    rmax = real(cd1.max) * z
    imin = imag(cd1.min) * z
    imax = imag(cd1.max) * z
    
    return ComplexDomain(rmin,rmax,imin,imax)
end;

function in(z::Complex,cd::ComplexDomain)
    return (real(cd.min) <= real(z) < real(cd.max) ) && (imag(cd.min) <= imag(z) < imag(cd.max) )
end;


##########################################







