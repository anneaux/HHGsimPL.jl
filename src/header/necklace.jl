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

### grads and hessians

	function grad(b::Beam, Ip::Float64,
	    q::Number,
	    ti::ComplexF64, tr::ComplexF64)

	    g = [dS_dti(b,Ip,q,ti,tr); dS_dtr(b,Ip,q,ti,tr)]
	    g = conj.(complex.(-1im .* g))
	    return g
	end

	function gradN(b::Beam, Ip::Float64,
	    q::Number,
	    ti::ComplexF64, tr::ComplexF64,
	    thresh::Float64 = 1.)

	    g = grad(b,Ip,q,ti,tr)
	    if norm(g) > thresh # bit lower than the gradient at the saddle point
	        return LinearAlgebra.normalize(g)
	    else 
	        return g
	    end
	end;

	function my_hessian(b::Beam, Ip::Float64,
	        q::Number,
	        ti::ComplexF64, tr::ComplexF64)

	    action(tvec) = real(-S_for_diff(b, Ip, tvec[1]+im*tvec[2],tvec[3]+im*tvec[4], q))
	    return FiniteDiff.finite_difference_hessian(action, [reim(ti)..., reim(tr)...])
	end


### point and lineseg
	mutable struct Point{T}
	    x::T
	    y::T
	    active::Bool
	end
	Point(x,y) = Point(x,y,true);


	mutable struct LineSeg
	    s::Union{UndefInitializer,Point{T}} where T<:Union{<:Number, Vector{<:Number}}
	    e::Union{UndefInitializer,Point{T}} where T<:Union{<:Number, Vector{<:Number}}
	    active::Bool
	    sindex::Union{UndefInitializer, Int}
	    eindex::Union{UndefInitializer, Int}
	end
	LineSeg(s::Point, e::Point) = LineSeg(s,e, s.active && e.active, undef,undef);
	LineSeg(s::Point, e::Point, active::Bool) = LineSeg(s, e, active, undef,undef);
	LineSeg(sindex::Int, eindex::Int, active::Bool) = LineSeg(undef, undef, active, sindex, eindex);
	LineSeg(sindex::Int, eindex::Int) = LineSeg(undef, undef, false, sindex, eindex);

    import Base.imag, Base.real
    imag(p::Point) = Point(imag.(p.x),imag.(p.y), p.active)
    real(p::Point) = Point(real.(p.x),real.(p.y), p.active);

    imag(ls::LineSeg) = LineSeg(imag(ls.s),imag(ls.e), ls.active)
    real(ls::LineSeg) = LineSeg(real(ls.s),real(ls.e), ls.active);

    function get_point(ls::LineSeg,which::Symbol=:s)
        if which==:s
            return Point(ls.s.x, ls.s.y)
        elseif which == :e
            return Point(ls.e.x, ls.e.y)
        else
            return println("You've got a problem!")
        end
    end



### necklacy things
function initialise!(necklace::Vector{LineSeg},points::Vector{Point},
        b::Beam, Ip::Float64,
        q::Number,
        ti::ComplexF64, tr::ComplexF64;
        Ninit::Int64 = 20,
        ϵ::Float64 = 0.01)

    hessian = my_hessian(b,Ip,q,ti,tr)

    # this could certainly be made more julian    
    eigenvectors = [[complex(vec[1:2]...), complex(vec[3:4]...)] for vec in eachcol(eigvecs(hessian))] 

    pointsini = ([[ti,tr] .+ ϵ * (cos(θ) * eigenvectors[3] + sin(θ) * eigenvectors[4]) for θ in range(0, stop=2π, length=Ninit+1)])

    push!(points, [Point(p[1], p[2]) for p in pointsini]...)
    push!(necklace, [LineSeg(i, i+1 , true) for i in 1:(length(points)-1)]...)
    push!(necklace, LineSeg(length(points), 1 , true)) # closing the necklace
end

### TODO this Δ could definitely get a more sophisticated default value
function subdivide!(lineseg::LineSeg,
            necklace::Vector{LineSeg}, points::Vector{Point};
            Δ::Float64=1.)    

    p1 = points[lineseg.sindex]
    p2 = points[lineseg.eindex]

    active = p1.active && p2.active

    Δx(p1::Point,p2::Point) = norm(p2.x - p1.x)
    Δy(p1::Point,p2::Point) = norm(p2.y - p1.y)
    midx(p1::Point,p2::Point) = (p2.x + p1.x)./2
    midy(p1::Point,p2::Point) = (p2.y + p1.y)./2

    if active && ( max(Δx(p1,p2), Δy(p1,p2)) > Δ)    
        lineseg.active = false # lineseg gets turned inactive when being divided.       
        midpoint = Point(midx(p1,p2), midy(p1,p2)) 
        push!(points, midpoint)

        lineseg1mid = LineSeg(lineseg.sindex,length(points),true)
        linesegmid2 = LineSeg(length(points),lineseg.eindex,true)
        push!(necklace, lineseg1mid)
        push!(necklace, linesegmid2)
    end
end

function flow!(necklace::Vector{LineSeg}, points::Vector{Point},
        b::Beam, Ip::Float64,
        q::Number;
        δ::Float64=0.1,
        threshold::Float64=0.5
        )

    for i in 1:length(points)
        # TODO check both real and imaginary part?
        if points[i].active # for the active points
            # set them to be active (= still flowing) if they are above threshold
            points[i].active = real(-im * S(b, Ip, points[i].x, points[i].y, q)) < 0 #(in Job's code that's h-function > thresh, I should clearly state which sign I'm using where etc.) 
            if points[i].active
                step = δ .* gradN(b, Ip, q, points[i].x, points[i].y, threshold)
                points[i].x += step[1]
                points[i].y += step[2]
            end
        end
    end
end

function adorn_necklace!(necklace::Vector{LineSeg}, points::Vector{Point})
    for i in 1:length(necklace)
        necklace[i] = LineSeg(points[necklace[i].sindex], points[necklace[i].eindex])
    end
end;

### get necklace
function get_necklace(b::Beam, Ip::Float64,
        q::Number,
        ti::ComplexF64, tr::ComplexF64
        ; Ninit::Int64=20, Ncounter::Int64=500,
        eigvecfactorinit::Float64 = 0.01, # I should come up with sophisticated guesses here.
        flowstepfactor::Float64 = 0.1, 
        subdividethreshold::Float64 = 0.5 )
       
    necklace = Vector{LineSeg}()
    points = Vector{Point}()

    initialise!(necklace, points, b, Ip, q, ti, tr, Ninit = Ninit, ϵ = eigvecfactorinit)

    ### find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(grad(b,Ip,q, p.x, p.y)) for p in points]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    
    counter = 0
        
    while counter < Ncounter
        counter += 1

        tmp = deepcopy(necklace)
        flow!(necklace, points, b, Ip, q, threshold = threshold, δ = flowstepfactor)

        if count([p.active for p in points]) == 0
            println("I broke because the flow stopped after $counter iterations")
            break
        end
        for i in 1:length(necklace)
            subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
        end
        keepat!(necklace, [ls.active for ls in necklace])
    end
        
    if counter == Ncounter 
        println("I broke because the counter reached its max, i.e. $Ncounter")
    end        

    adorn_necklace!(necklace, points)
    
    return necklace
end;
