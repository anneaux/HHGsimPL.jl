### everything related to saddle points

### TODO implement different methods to make sure all SPs are found: classical guess, sobol, nroots



abstract type SP end

### Saddle ###########################
struct Saddle <: SP
	q::T where T <: Number
	ti::Complex{Float64}
	tr::Complex{Float64}
	p::Vector{ComplexF64} # do I want this here???
end

Base.show(io::IO, s::Saddle) = print(io,
    "ti: $(round(s.ti,sigdigits=5)), tr: $(round(s.tr,sigdigits=5)), q$(round(s.q,sigdigits=5))")

import Base./
function /(s::Saddle, TC::Real)
    return Saddle(s.q, s.ti/TC, s.tr/TC, s.p)
end


#########################################################

function solve_SPEqs(q::Number, t0::Vector{T}, b::Beam, Ip::Real,
    roundDigits::Int64=5) where T <: Real

    try
    ### using NLsolve
        function speqs!(F, x)
            F[1] = real(speq1(b, Ip, x[1], x[2], x[3], x[4]))
            F[2] = real(speq2(b, Ip, q, x[1], x[2], x[3], x[4]))
            F[3] = imag(speq1(b, Ip, x[1], x[2], x[3], x[4]))
            F[4] = imag(speq2(b, Ip, q, x[1], x[2], x[3], x[4]))
        end
        
        result = nlsolve(speqs!, t0)#, method = :trust_region, factor =fac )#, ftol = 1e-13)#, method = :anderson)
        if converged(result)
            tiSP = result.zero[1] + im*result.zero[2]
            trSP = result.zero[3] + im*result.zero[4]
            tiSP = round(tiSP, digits = roundDigits)
            trSP = round(trSP, digits = roundDigits)
            return tiSP, trSP
        else
            return nothing, nothing
        end
    catch e
        println("Error in solve_SPEqs(): $e")
        return nothing,nothing
    end
end


function check_sp(b::Beam, 
    ti::Complex{Float64}, tr::Complex{Float64};
    tt_minimal::Float64=0.005)

    tt = tr - ti # travel time

	# this is just the 'boundary condition' that recombination happens after ioniosation (i.e., tr > ti)
    check1 = tt_minimal*TCycle(b) <  real(tt)

    # and we only look at positive tunnelling times
    check2 = imag(ti) > 0 
    return check1 && check2
end


function check_sp(b::Beam, 
    ti::Nothing, tr::Nothing;
    tt_minimal::Float64=0.005)

    return false
end


#########################################################

function find_saddles_sobol(q::Number; 
        b::Beam, Ip::Real,
        ti_cd::ComplexDomain, 
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min) - imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        N::Int64=200, # number of seeds generated per domain
        tt_minimal::Float64 = 0.005,
        beam::Beam=b ) 
    roundDigits = 2 # I should certainly think this over it seems too much

	saddles = Vector{Saddle}()
    
    ti_seq = SobolSeq(reim(ti_cd.min),reim(ti_cd.max))
    tr_seq = SobolSeq(reim(tr_cd.min),reim(tr_cd.max))

	for i in 1:N
        ti0 = next!(ti_seq)
        tr0 = next!(tr_seq)
        
        t0 = [ti0[1]; ti0[2]; tr0[1]; tr0[2]]

        tiSP, trSP = solve_SPEqs(q, t0, beam, Ip, roundDigits) 

        ### check conditions and deposit in array
        if check_sp(b, tiSP,trSP, tt_minimal = tt_minimal) == true && 
            in(tiSP, ti_cd) && # maybe I want this to be an option
            in(trSP,tr_cd)
            # &&new
            push!(saddles, Saddle(q, tiSP, trSP, p_stationary(b, tiSP, trSP)))
        end
	end

    unique!( s -> round.([s.tr,s.ti], digits=roundDigits), saddles)

    sort!(saddles, by = x -> real(x.ti))
	return saddles
end;


# I don't actually need all those kwargs but I want to make it easy to switch between the functions and I don't know how to do that nicely for now
function find_saddle_similar_seed(q::Number, seed::Saddle;
        b::Beam, Ip::Real,
        ti_cd::ComplexDomain,
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min) - imag(ti_cd.max)*im,ti_cd.max + TCycle(b)), 
        tt_minimal::Float64 = 0.005,
        beam::Beam=b ) # number of seeds generated per domain
    
    roundDigits = 2
    t0 = [reim(seed.ti)... , reim(seed.tr)...]

    tiSP, trSP = solve_SPEqs(q, t0, beam, Ip, roundDigits) # last arg: rounddigits

    if check_sp(b, tiSP, trSP, tt_minimal = tt_minimal) && in(tiSP, ti_cd) && in(trSP,tr_cd)
        return Saddle(q, tiSP, trSP, p_fun(b, tiSP,trSP))
    else 
        return nothing
    end
end;

# I could then somehow include it in a function like this
#         previous_saddles = saddles_Dict[q]
#         next_saddles = Vector{Saddle}()
#         q += 1
# #         @show q
#         for s in previous_saddles
#             new = find_saddle_similar_seed(q, s, 
#                 b = beam, Ip = Ip, ti_cd = ti_domain, tr_cd = tr_domain, N = 200)   
#             if new != nothing
#                 push!(next_saddles, new)
#             else 
#                 println("couldn't find a new one!")
#             end
#         end
#         unique!( s -> round.([s.tr,s.ti], digits= 2), next_saddles)
#         sort!(next_saddles, by = x -> real(x.ti))
#         saddles_Dict[q] = next_saddles