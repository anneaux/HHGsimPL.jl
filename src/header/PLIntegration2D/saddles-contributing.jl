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


    function find_crossing(line::Vector{LineSeg}, point::Point{T}, tolerance::Float64=1.;
        loginfo=[]) where T<:Real

        distances = [distance_point_to_line(point, seg) for seg in line]
        
        # finds local minima of the distances, filters for those where the height is <0.8, and returns the respective indices
        # https://docs.juliahub.com/Peaks/3TWUM/0.5.2/
        intersections = findminima(vcat(distances, distances[1:min(20, length(distances))])) |> peakheights(;max = tolerance) |> peakproms(;min = 0.5)
        peakindices =  unique(mod1.(intersections.indices, length(distances)))

        ### double-check that peaks are smaller than norm, think: adaptive tolerance for peak height. averaging over norms in that region because otherwise sometimes I'm unlucky

        filter!(pidx -> distances[pidx] < sum([norm(ls) for ls in line[mod1.(collect(pidx-2:pidx+2), length(line))]])/5, peakindices)
 
        if length(peakindices) == 1
           return peakindices[1]
        elseif length(peakindices) == 0
            return nothing
        else
            @warn "I'm hitting the integration plane more than once I think"
            # log_error("new-necklace-hitting-ID-errors.txt", "Warning (2) for $(loginfo).")
            return peakindices[1]
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
function real_projected_contourlines(
    f::Function,
    ti::ComplexF64, tr::ComplexF64,
    ti_range::Real=50, tr_range::Real=50
    ; Ntimes = 100)    
    
    tir_values = range(real(ti)- ti_range, stop = real(ti) + ti_range, length = Ntimes)
    trr_values = range(real(tr)- tr_range, stop = real(tr) + tr_range, length = Ntimes)
 
    ### level line for the saddle point
    S_values = [f(complex(tir), complex(trr)) for tir in tir_values, trr in trr_values]
    S_saddle = f(ti, tr)
    contour_saddle = Contour.contour(tir_values, trr_values, imag.(S_values), imag(S_saddle) )

    return contour_saddle.lines
end

### checking if conditions are fulfilled
function check_contribution(necklace::Vector{LineSeg}, 
    f::Function,
    ti::ComplexF64, tr::ComplexF64,
    ti_range::Real=50, tr_range::Real=50
    ; Ntimes = 100 )
    
    ### check if necklace hits real plane
    p = Point(0.,0.)
    idx = find_crossing( imag.(necklace), p) # can add loginfo here
    
    if isnothing(idx)
        @debug "it doesn't contribute! (1)"
       active = false
    else         
        ### get the point where it hits & check if it's in the integration domain
        hitting_point = real(necklace[idx].s.y) > real(necklace[idx].s.x) ? 
            get_point(real(necklace[idx])) : nothing
        # mustn't use the starting point here, could use the centre point!
        
        if isnothing(hitting_point)
           println("it doesn't contribute! (2)") # because this shouldn't happen!
           active = false
        else
            ### check if the projected contour runs through that point
            contourlines = real_projected_contourlines(f, ti, tr, ti_range, tr_range)
            
            crosses = [false]
            for line in contourlines
                crossings = find_crossing(line, hitting_point)
                @debug "crosses at $crossings"
                push!(crosses, !isnothing(crossings))
            end
            active = any(crosses)
            if !active @debug "it doesn't contribute! (3)" end
        end
    end       

    return active
end;

function check_contribution(necklace::Nothing, 
    f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64,
    ti_range::Real=50, tr_range::Real=50
    ; Ntimes = 100 )
    return false
end

function check_contribution(necklace::Nothing, 
    f::Function,
    ti::ComplexF64, tr::ComplexF64,
    ti_range::Real=50, tr_range::Real=50
    ; Ntimes = 100 )
    return false
end


function check_contribution(
    f::Function,
    f_grad::Function,
    f_hessian::Function,
	ti::ComplexF64, tr::ComplexF64,
    ti_range::Real=50, tr_range::Real=50
    ; Ntimes::Int64 = 100, logerrors::Bool=false, kwargs...)
    
    if real(f(ti, tr)) < 0
        necklace = get_necklace(f,f_grad,f_hessian, ti, tr; logerrors=logerrors, kwargs...)
        check_contribution(necklace, f, ti, tr, ti_range, tr_range, Ntimes = Ntimes)
    else 
        @debug "it doesn't contribute! (0)"
        return false
    end
    # what happens if doesn't converge?   
end

nothing