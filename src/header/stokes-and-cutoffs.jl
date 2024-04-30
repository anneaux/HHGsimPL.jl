### everything related to the cutoff and the Stokes phenomena

# TODO I should make this more general, because not every fold point should be referred to as a 'harmonic cutoff'. E.g what about the 'lower cutoffs'?

###### Thc ########### 
# harmonic cutoff times, i.e. solutions to the second derivative
struct Thc
  tihc::ComplexF64
  trhc::ComplexF64
  qhc::ComplexF64

  #Thc() = new()
  #Thc(vec::Vector{ComplexF64}) = new(vec[1],vec[2])
end


Base.show(io::IO, thc::Thc) = println(io,
    "tihc: $(round((thc.tihc),sigdigits=5)), trhc: $(round((thc.trhc),sigdigits=5)), qhc: $(round((thc.qhc),sigdigits=5))")


### thc equations ###########

function thceq1(b::Beam, Ip::Float64,
  ti::ComplexF64, tr::ComplexF64) ### derived from eq. 15 in Emilio's paper

  return d2Sv_dtr2(b, Ip, ti,tr) * d2Sv_dti2(b, Ip, ti,tr) - d2Sv_dtitr(b, Ip, ti,tr) * d2Sv_dtitr(b, Ip, ti,tr)
end




function thceq2(b::Beam, Ip::Float64,
  ti::ComplexF64,tr::ComplexF64)

  return dSv_dti(b, Ip, ti,tr)
end   


### cutoff energy 
function E_hc(b::Beam, Ip::Float64,
  tihc::ComplexF64,trhc::ComplexF64)
  return dSv_dtr(b, Ip, tihc, trhc)
end
E_hc(b::Beam, Ip::Float64, thc::Thc) = E_hc(b,Ip,thc,tihc,thc.trhc)

### harmonic order 
function thc_qc(b::Beam, Ip::Float64,
  tihc::ComplexF64,trhc::ComplexF64)
  return E_hc(b,Ip,tihc,trhc)/b.omega1
end
thc_qc(b::Beam, Ip::Float64, thc::Thc) = thc_qc(b, Ip, thc.tihc, thc.trhc)




#########################################################

function solve_Thceqs(t0::Vector{T}, 
    b::Beam, Ip::Real,
    roundDigits::Int64=5) where T <: Real

    try
    ### using NLSolve
        #
        function thceqs!(F,x)

            F[1] = real(thceq1(b, Ip, x[1]+x[2]*im, x[3]+x[4]*im))
            F[2] = real(thceq2(b, Ip, x[1]+x[2]*im, x[3]+x[4]*im))
            F[3] = imag(thceq1(b, Ip, x[1]+x[2]*im, x[3]+x[4]*im))
            F[4] = imag(thceq2(b, Ip, x[1]+x[2]*im, x[3]+x[4]*im))

        end

        result = nlsolve(thceqs!, t0)#, method = :trust_region, factor =fac )#, ftol = 1e-13)#, method = :anderson)
        if converged(result)
            tiSP = result.zero[1] + im*result.zero[2]
            trSP = result.zero[3] + im*result.zero[4]
            tiSP = round(tiSP, digits = roundDigits)
            trSP = round(trSP, digits = roundDigits)
            return tiSP, trSP
        else 
            return nothing,nothing
        end
    catch e
        println("Error in solve_Thceqs(): $e")
        return nothing,nothing
    end
end


function find_thcs_sobol(b::Beam, Ip::Real,
        ti_cd::ComplexDomain, 
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min) - imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        N::Int64=20; tt_minimal::Float64 = 0.15,
        beam::Beam=b # just for clarity
        )
    
    thc_array = Vector{Thc}()

    ti_seq = SobolSeq(reim(ti_cd.min),reim(ti_cd.max))
    tr_seq = SobolSeq(reim(tr_cd.min),reim(tr_cd.max))
    
    for i in 1:N
        tihc0 = next!(ti_seq)
        trhc0 = next!(tr_seq)
        
    thc0 = [tihc0[1]; tihc0[2]; trhc0[1]; trhc0[2]]

        tihc, trhc = solve_Thceqs(thc0, beam, Ip, 5)

        ### checking if already found this
        new = length(findall(x -> x.tihc == tihc && x.trhc == trhc, thc_array)) == 0     # unique
        if check_sp(b, tihc,trhc,tt_minimal = tt_minimal) && 
            new && in(tihc,ti_cd) && in(trhc,tr_cd)
            qhc = thc_qc(b, Ip, tihc, trhc)
            push!(thc_array, Thc(tihc, trhc, qhc))       
        end 
    end 
    
    sort!(thc_array, by= x -> real(x.tihc))
    return thc_array
end;
