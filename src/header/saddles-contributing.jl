### everything to decide whether or not a given saddle point contributes.
### this could be implemented in various methods again. Also maybe it should give a warning if there're multiple saddle points nearby and if a Gaussian approximation is a bad idea?




### utils for deciding whether a line crosses a given point
	function distance_point_to_line(p::AbstractVector, s::AbstractVector, t::AbstractVector)
	    midpoint = (s .+ t)./2
	    return norm(p .- midpoint)
	end

	function distance_point_to_line(p::Point, l::LineSeg)
	    return distance_point_to_line([p.x,p.y], [l.s.x, l.s.y], [l.e.x, l.e.y])
	end

	function find_crossing(line::Vector{LineSeg}, point::Point{T}, tolerance::Float64=0.8) where T<:Real
	    mindist, index = findmin([distance_point_to_line(point, seg) for seg in line])
	    if mindist < tolerance
	        return index
	    else 
	        return nothing
	    end
	end

	function find_crossing(curve::Curve2{Tuple{T, T}}, point::Point{T}, tolerance::Float64=0.8) where T<:Real
	    line = [LineSeg( Point(curve.vertices[i]...), Point(curve.vertices[i+1]...)) for i in 1:(length(curve.vertices)-1) ]
	    return find_crossing(line, point, tolerance)
	end

	function find_crossing(nocurve::Missing, point::Point{T}, tolerance::Float64=0.8) where T<:Real
	    return nothing
	end

### calculating the contour line through a given saddle
function contourline_through_saddle(b::Beam, Ip::Float64,
    q::Number,
    ti::ComplexF64, tr::ComplexF64,
    ti_cd::ComplexDomain, tr_cd::ComplexDomain
    ; Ntimes = 100)    
    
    TC = TCycle(b)
    tir_values = range(real(ti)- 0.5TC, stop = real(ti) + 0.5TC, length = Ntimes)
    tii_values = range(-1., stop = imag(ti) + 0.25TC, length = Ntimes)
    trr_values = range(real(tr)- 0.5TC, stop = real(tr) + 0.5TC, length = Ntimes)
    tri_values = range(-1., stop = imag(ti) + 0.25TC, length = Ntimes) 

    ### level line for the saddle point
    S_values = [-1im*S(b, Ip, complex(tir), complex(trr), q) for tir in tir_values, trr in trr_values]
    S_saddle = -1im*S(b, Ip, ti, tr, q)
    contour_saddle = Contour.contour(tir_values, trr_values, imag.(S_values), imag(S_saddle) )

    if length(contour_saddle.lines) == 1
        return contour_saddle.lines[1]
    elseif length(contour_saddle.lines) == 0
        println("Careful! There's no level line going through the saddle point for $b at q $q.")
        return missing
    else
        print("Careful! There's more than one level line going through the saddle point for $b at q $q,")
        p = Point(real(ti), real(tr))
        filter!(curve -> (find_crossing(curve, p, 4.) != nothing ), contour_saddle.lines) # where does this value come from and should I choose a better one?
        if length(contour_saddle.lines) == 1
            println(" but I've resolved it.")   
            return contour_saddle.lines[1]
        else 
            println(" and we do in fact have a problem.") #. There was more than one level line going through the saddle point for $b at q $q, but they seem to not actually intersect (or multiple of them)")
        end
    end
end

function contourline_through_saddle(b::Beam, Ip::Float64,
	s::Saddle,
    ti_cd::ComplexDomain, tr_cd::ComplexDomain
    ; Ntimes = 100) 

    contourline_through_saddle(b, Ip, s.q, s.ti, s.tr, ti_cd, tr_cd; Ntimes = Ntimes) 
end

### checking if conditions are fulfilled
function check_contribution(necklace::Vector{LineSeg}, 
    b::Beam, Ip::Float64,
    q::Number,
    ti::ComplexF64, tr::ComplexF64,
    ti_cd::ComplexDomain, tr_cd::ComplexDomain
    ; Ntimes = 100 )
    
    ### check if necklace hits real plane
    p = Point(0.,0.)
    idx = find_crossing( imag.(necklace), p)
    
    if isnothing(idx)
       println("it doesn't contribute! (1)")
       active = false
    else         
        ### get the point where it hits & check if it's in the integration domain
        hitting_point = real(necklace[idx].s.y) > real(necklace[idx].s.x) ? 
            get_point(real(necklace[idx])) : nothing
        # mustn't use the starting point here, could use the centre point!
        
        if isnothing(hitting_point)
           println("it doesn't contribute! (2)")
           active = false
        else
            ### check if the contour runs through that point
            contourline = contourline_through_saddle(b, Ip, q, ti, tr, ti_cd, tr_cd )
            
            if isnothing(find_crossing(contourline, hitting_point))
                println("it doesn't contribute! (3)")
                active = false
            else 
                println("it contributes!")
                active = true
            end
        end
    end
    return active
end


function check_contribution(b::Beam, Ip::Float64,
	q::Number,
	ti::ComplexF64, tr::ComplexF64,
    ti_cd::ComplexDomain, tr_cd::ComplexDomain
    ; Ntimes = 100, Ncounter = 600)
    
    if real(-im * S(b, Ip, ti, tr, q)) < 0
        necklace = get_necklace(b, Ip, q, ti, tr, Ncounter = Ncounter)
        check_contribution(necklace, b, Ip, q, ti, tr, ti_cd, tr_cd, Ntimes = Ntimes)
    else 
        return false
    end
    # what happens if doesn't converge?   
end