### Constants and Units

const EAU = 27.2114
const IAU = 3.5*10^16
const LAU = 0.052918
const TAU = 2.419*10^-17
const alpha = 1. /137
const c = 1/alpha

const cNU = 299792458;


natural_time(atomictime::Float64) = atomictime * LAU * c  * 1e-9 /cNU
to_fs(atomic_time::Float64) = natural_time(atomic_time)*1e15;

IpAU(IpSI::Real) = IpSI/EAU

####
