
# using Distances # das brauch ich!
using DataFrames

include("meanshift-Anne.jl")
using .MeanShiftClustering

# requires HHGsimPL.jl but I don't know how to do that



function cluster_ionisation_times(vals::Vector{T}) where T<:Real
#     vals = test_x
    msr = MeanShiftClustering.meanshift(vals)
    
    ionisationbursts_df = DataFrame(index = Int64[],
        center=Float64[],
        range = Tuple{Number,Number}[],
        data = Vector{Number}[])
    
    cent_prev = minimum(vals)

    if (msr.centers[1,:]) == sort(msr.centers[1,:])
        centers = msr.centers[1,:]
        groupindices = msr.assignments
    else 
        # println("had to sort the centers")
        sorted_indices = sortperm(msr.centers[1,:])

        # Create a mapping from old indices to new indices
        index_mapping = Dict(sorted_indices[i] => i for i in 1:length(sorted_indices))

        # Sort the centers
        centers = sort(msr.centers[1,:])

        # Update the assignments using the index mapping
        groupindices = [index_mapping[a] for a in msr.assignments]        
    end

    for (ci, c) in enumerate(centers)
        @assert (c>cent_prev) "the ionisationbursts are not in order!"
        
        cent_prev = c
        data = vals[groupindices .== ci]
        push!(ionisationbursts_df, [ci, c, (minimum(data), maximum(data)), data ])
    end
    return ionisationbursts_df
end



function get_saddles_dict(b::Beam, Ip::Real, q_values::Vector{T}, 
    ti_cd::ComplexDomain,   
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min) - imag(ti_cd.max)*im,ti_cd.max + TCycle(b));
    Nsobolseeds::Int64 = 2000,
    filterfun::Function = s -> 10 < real(s.tr - s.ti) < TCycle(b),
    check_contrib::Bool = true) where T <: Real
    
    saddlesTup_Dict = Dict{Real,Vector{Tuple{Saddle,Bool}}}()
    
    # Threads.@threads 
    for q in q_values
#         saddles = find_saddles_sobol(q, b=b, Ip = Ip, ti_cd = ti_cd, tr_cd = tr_cd, N = 2000)
        saddles = filter(filterfun, find_saddles_sobol(q, b=b, Ip = Ip, ti_cd = ti_cd, tr_cd = tr_cd, N = Nsobolseeds) )
            
        vec = Vector{Tuple{Saddle,Bool}}()
        for s in saddles
            relevant = check_contrib ? check_contribution(b, Ip, q, s.ti, s.tr, ti_cd, tr_cd, Ncounter = 1000) : true 
            push!(vec, (s,relevant))
        end
        saddlesTup_Dict[q] = vec
    end 
    
    return saddlesTup_Dict
    
end;



function get_ionisationbursts(beam::Beam, Ip::Number, q_values::Vector{T}) where T <: Real
        
    ti_domain = ComplexDomain(0.,1., 0.01,0.35)*TCycle(beam);
    tr_domain = ComplexDomain(0. +0.1, 1. + 1., -0.35,0.35)*TCycle(beam);
    
    saddles_ini = get_saddles_dict(beam, Ip, q_values, ti_domain, tr_domain, Nsobolseeds = 1000
        , filterfun = s -> 10 < real(s.tr - s.ti) < 1*TCycle(beam), check_contrib = false)
    
    data = [real(sTup[1].ti)/TCycle(beam) for sTup in vcat(values(saddles_ini)...)]   
    
    ibs_df = cluster_ionisation_times(data)

    lims = Vector{Float64}()
    for row in eachrow(ibs_df)  
        next_min = if row.index !==(size(ibs_df)[1])
            ibs_df[row.index+1,:range][1]
            else
                ibs_df[1,:range][1] + 1.
            end 
        midpoint = (next_min + row.range[2])/2
        push!(lims, midpoint)
    end
    push!(lims, lims[end]-1.)
    sort!(lims)

    ibs_df[!, "limits"] = collect(zip(lims[1:end-1], lims[2:end]))

    return ibs_df
end

### getting ionisationburst limits quicker
ib_limits_dict = Dict{Beam, Vector{Tuple{Float64, Float64}}}()

function get_ionisationburst_limits(beam::Beam, Ip::Number, q_values::Vector{T}=collect(15:2:40.)) where T <: Real

    if haskey(ib_limits_dict, beam)
        return ib_limits_dict[beam]
    else
        try
            ionisationbursts = get_ionisationbursts(beam, Ip, q_values)
            lims = [r.limits for r in eachrow(ionisationbursts)]
            ib_limits_dict[beam] = lims
            return lims
        catch e
            println("Error in get_ionisationburst_limits(): ", e)
            lims = [(0., 0.5),(0.5, 1.0)]
            ib_limits_dict[beam] = lims
            return lims
        end
    end
end

nothing