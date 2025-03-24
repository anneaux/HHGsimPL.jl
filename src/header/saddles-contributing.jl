### everything to decide whether or not a given saddle point contributes.
### this could be implemented in various methods again. Also maybe it should give a warning if there're multiple saddle points nearby and if a Gaussian approximation is a bad idea?

dist(p1::Point, p2::Point) = norm(@SVector[p2.x-p1.x,p2.y-p1.y])
import Base.length 
length(ls::LineSeg) = dist(ls.e, ls.s)


### utils for deciding whether a line crosses a given point
	function distance_point_to_line(p::AbstractVector, s::AbstractVector, t::AbstractVector)
	    midpoint = (s .+ t)./2
	    return norm(p .- midpoint)
	end

	function distance_point_to_line(p::Point, l::LineSeg)
	    return distance_point_to_line([p.x,p.y], [l.s.x, l.s.y], [l.e.x, l.e.y])
	end

	# function find_crossing(line::Vector{LineSeg}, point::Point{T}, tolerance::Float64=0.8) where T<:Real
	#     mindist, index = findmin([distance_point_to_line(point, seg) for seg in line])

	#     if mindist < tolerance
	#         return index
	#     else 
	#         return nothing
	#     end
	# end


    function average_distance(line::Vector{LineSeg}, pidx::Int64, threshold::Float64=0.5) # flowstepfactor
        line_region = [line[pidx]]

        for r in 1:min(10, length(line)-pidx-1)
            push!(line_region, line[pidx+r])        
            if norm(line[pidx+r]) > threshold # 2*flowstepfactor
                break
            end
        end   
        for r in -1:-1:-min(10, pidx-1)
            push!(line_region, line[pidx+r])
            if norm(line[pidx+r]) > threshold
                break
            end
        end    
        
        av_dist = sum([norm(ls) for ls in line_region])/length(line_region)#nregion
    #     @show av_dist
        return av_dist
    end






    function find_crossing(line::Vector{LineSeg}, point::Point{T}, tolerance::Float64=1.; threshold::Float64=0.5,
        loginfo=[]) where T<:Real

        distances = [distance_point_to_line(point, seg) for seg in line]
        
        # finds local minima of the distances, filters for those where the height is < tolerance, and returns the respective indices
        # https://docs.juliahub.com/Peaks/3TWUM/0.5.2/
        intersections = findminima(vcat(distances, distances[1:min(20, length(distances))])) |> peakheights(;max = tolerance) |> peakproms(;min = 0.5)
        peakindices =  unique(mod1.(intersections.indices, length(distances)))

        ### double-check that peaks are smaller than norm, think: adaptive tolerance for peak height. averaging over norms in that region because otherwise sometimes I'm unlucky
        # region = 3
        # nregion = 2*region + 1

        # filter!(pidx -> distances[pidx] < sum([norm(ls) for ls in line[mod1.(collect(pidx-region:pidx+region), length(line))]])/nregion, peakindices)
 
        filter!(pidx -> distances[pidx] < average_distance(line, pidx, threshold), peakindices)

        if length(peakindices) == 1
           return peakindices[1]
        elseif length(peakindices) == 0
            return nothing
        else
            @warn "I'm hitting the integration plane more than once I think"
            log_error("new-necklace-hitting-ID-errors.txt", "Warning (2) for beam $(loginfo[1]) at q $(loginfo[2]) with ti $(loginfo[3]) and tr $(loginfo[4]).")
            return peakindices[1]
        end
    end


	function find_crossing(curve::Curve2{Tuple{T, T}}, point::Point{T}, tolerance::Float64=0.8; threshold::Float64=0.5) where T<:Real
	    line = [LineSeg( Point(curve.vertices[i]...), Point(curve.vertices[i+1]...)) for i in 1:(length(curve.vertices)-1) ]
	    return find_crossing(line, point, tolerance, threshold=threshold)
	end

	function find_crossing(nocurve::Missing, point::Point{T}, tolerance::Float64=0.8; threshold::Float64=0.5) where T<:Real
	    return nothing
	end

### calculating the contour line through a given saddle
# function real_projected_contourlines(b::Beam, Ip::Float64,
#     q::Number,
#     ti::ComplexF64, tr::ComplexF64,
#     ti_cd::ComplexDomain, tr_cd::ComplexDomain
#     ; Ntimes = 100)    
    
#     TC = TCycle(b)
#     tir_values = range(real(ti)- 0.5TC, stop = real(ti) + 0.5TC, length = Ntimes)
#     # tii_values = range(-1., stop = imag(ti) + 0.25TC, length = Ntimes)
#     trr_values = range(real(tr)- 0.5TC, stop = real(tr) + 0.5TC, length = Ntimes)
#     # tri_values = range(-1., stop = imag(ti) + 0.25TC, length = Ntimes)  # this is wrong, because tr can have negative imaginary part! Luckily I don't need that here anyway ;-)

#     ### level line for the saddle point
#     S_values = [-1im*S(b, Ip, complex(tir), complex(trr), q) for tir in tir_values, trr in trr_values]
#     S_saddle = -1im*S(b, Ip, ti, tr, q)
#     contour_saddle = Contour.contour(tir_values, trr_values, imag.(S_values), imag(S_saddle) )

#     return contour_saddle.lines
# end

# function real_projected_contourlines(b::Beam, Ip::Float64,
#     s::Saddle,
#     ti_cd::ComplexDomain, tr_cd::ComplexDomain
#     ; Ntimes = 100) 

#     real_projected_contourlines(b, Ip, s.q, s.ti, s.tr, ti_cd, tr_cd; Ntimes = Ntimes) 
# end

### checking if conditions are fulfilled
function check_contribution(necklace::Vector{LineSeg}, 
    b::Beam, Ip::Float64,
    q::Number,
    ti::ComplexF64, tr::ComplexF64
    ; Ntimes = 100, kwargs...)
    # flowstepfactor::Float64=0.5 )
    
    flowstepfactor = try kwargs[:flowstepfactor] catch e 0.8 end
    ### check if necklace hits real plane
    p = Point(0.,0.)
    idx = find_crossing( imag.(necklace), p, loginfo=(b,q,ti,tr), threshold = 2*flowstepfactor)
    
    if isnothing(idx)
        @debug "it doesn't contribute! (1)"
       active = false
    else         
        ### get the point where it hits & check if it's in the integration domain
        hitting_point = real(necklace[idx].s.y) > real(necklace[idx].s.x) ? 
            get_point(real(necklace[idx])) : nothing
        # mustn't use the starting point here, could use the centre point!
        
        if isnothing(hitting_point)
           println("it doesn't contribute! $q (2)") # because this shouldn't happen!
           active = false
        else
            H_at_hp = imag(-1im*S(b, Ip, necklace[idx].s.x, necklace[idx].s.y, q ))
            H_at_sp = imag(-1im*S(b, Ip, ti,tr, q))
            if abs(H_at_hp - H_at_sp) < 1.
                active = true
            else
                println("it doesn't contribute! (4)")
                active = false
            end
        end
    end       

    return active
end;

function check_contribution(necklace::Nothing, 
    b::Beam, Ip::Float64,
    q::Number,
    ti::ComplexF64, tr::ComplexF64
    ; Ntimes = 100, kwargs... )
    return false
end



function check_contribution(b::Beam, Ip::Float64,
	q::Number,
	ti::ComplexF64, tr::ComplexF64
    ; Ntimes::Int64 = 100, logerrors::Bool=false, kwargs...)
    # Ncounter = 600, logerrors::Bool=false)
    
    if real(-im * S(b, Ip, ti, tr, q)) < 0
        necklace = get_necklace(b, Ip, q, ti, tr; logerrors=logerrors, kwargs...)
        check_contribution(necklace, b, Ip, q, ti, tr; Ntimes = Ntimes, kwargs...)
    else 
        @debug "it doesn't contribute! (0)"
        return false
    end
    # what happens if doesn't converge?   
end