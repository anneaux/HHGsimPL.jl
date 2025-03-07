### relevant functions to determine crossing etc.
struct Line0
    s::StaticArrays.SVector{2,Float64}
    e::StaticArrays.SVector{2,Float64}

    Line0(t1::Tuple{Float64, Float64}, t2::Tuple{Float64, Float64}) =
        new(SA[t1...],SA[t2...])
end

function intersection(l1::Line0, l2::Line0) 
    x1, y1 = l1.s
    x2, y2 = l1.e
    x3, y3 = l2.s
    x4, y4 = l2.e
    
    denom = (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1)
    
    if denom == 0
        return nothing  # Line0s are parallel
    end
    
    ua = ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom

    ub = ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom
    
    if 0 <= ua <= 1 && 0 <= ub <= 1
        intersection = (x1 + ua*(x2 - x1), y1 + ua*(y2 - y1))
        return Point(intersection...)  # Line0s intersect
    else
        return nothing  # Line0s do not intersect
    end
end

function intersection(c1::Curve2,c2::Curve2)
    intersection_points = Vector{Point}()
    for i in 1:(length(c1.vertices)-1), i2 in 1:(length(c2.vertices)-1)
        l1 = Line0(c1.vertices[i],c1.vertices[i+1])
        l2 = Line0(c2.vertices[i2],c2.vertices[i2+1])
        p = intersection(l1,l2)
        if p != nothing
            push!(intersection_points,p)
        end 
         
    end 
    return intersection_points
end 

function intersection(con1_curves::ContourLevel, con2_curves::ContourLevel)
    intersection_points = Vector{Point}()
    for curve1 in lines(con1_curves)
        for curve2 in lines(con2_curves)
            push!(intersection_points,intersection(curve1,curve2)...)
        end 
    end 
    return intersection_points
end



function crosses_point(line::Line0, point::Point{T}, tolerance::Float64=0.25) where T<:Real
    if line.s.x <= line.e.x
        sx = line.s.x
        ex = line.e.x
    else
        ex = line.s.x
        sx = line.e.x
    end

    if line.s.y <= line.e.y
        sy = line.s.y
        ey = line.e.y
    else
        ey = line.s.y
        sy = line.e.y
    end

    sx -= tolerance*abs(sx-ex)
    ex += tolerance*abs(sx-ex)
    sy -= tolerance*abs(sy-ey)
    ey += tolerance*abs(sy-ey)

    xbound = (sx <= point.x <= ex)
    ybound = (sy <= point.y <= ey)
    return (xbound && ybound)
end


function crosses_point(c::Curve2, p::Point{Float64}, tolerance::Float64=0.25)
    crossing = false
    for i in 1:(length(c.vertices)-1)
        l = Line0(c.vertices[i],c.vertices[i+1])
        if crosses_point(l, p , tolerance)
            crossing = true
        end
    end
    return crossing
end

function is_saddle_relevant_ti(b::Beam, Ip::Real, saddle::Saddle,
    ti_cd::ComplexDomain = ComplexDomain(),
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b))
    )
    
    tr = saddle.tr 
    TC = TCycle(b)
    
    tir_values = range(real(ti_cd.min)-0.1*TC, stop=real(ti_cd.max)+0.1*TC, length=51) # I'm adding the ±0.1 because the contour line might reach out tout the domain otherwise
    tii_values = range(
        min(imag(ti_cd.min), 0.), 
        stop = imag(ti_cd.max), 
        length=51)    
    
    Sv_values = [S_v(b, Ip, ti, tr) for ti in (tir_values' .+ im*tii_values)]
  
    Sv_saddle = S_v(b, Ip, saddle)
    relevant = false
    
    mycurve = Curve2([
            (minimum(tir_values)./TC, 0.),
            (maximum(tir_values)./TC, 0.)])   # this is my original integration contour, agaiit's relevant I account for ausschweifende contour lines
    
    # original integration domain (line)
    intdomain_ti = Curve2([(minimum(tir_values)./TC, 0.), (maximum(tir_values)./TC, 0.)])    
  
    con_S_saddle_real = Contour.contour(tir_values./TC, tii_values./TC, real.(Sv_values'), real.(Sv_saddle))
    
    #     println(imag(Sv_saddle))
    for curve in lines(con_S_saddle_real)
        ip = intersection(curve, intdomain_ti)
        for ipx in ip
            ti_ip = (ipx.x + ipx.y*im)*TC
            Sv_ip = S_v(b, Ip, ti_ip, tr)
    #             println(imag(Sv_ip))
            ascending = round(imag(Sv_ip - Sv_saddle),digits=5) >= 0
            crossing = crosses_point(curve,Point(reim(saddle.ti)./TC...),max(tii_values[2]-tii_values[1],1.5))
            if ascending && crossing
                relevant = true
    #                 println("yes, there's a good ti one")
                return relevant
            end
        end
        
    end 
    return relevant
end;


function is_saddle_relevant_tr(b::Beam, Ip::Real, saddle::Saddle,
    ti_cd::ComplexDomain = ComplexDomain(),
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b))
    )
    
    TC = TCycle(b)
    if abs(imag(saddle.tr)) <= 0.01*TC
        # println("I'm in regime (1), because ", abs(imag(saddle.tr)) )
        return true
        
    elseif abs(imag(saddle.tr)) < 0.05*TC # if the SP is nearby the real axis
        # println("I'm in regime (2), because ", abs(imag(saddle.tr)) )
        str = saddle.tr
        # when the sp is exactly at the border I want to make sure it reaches up to line
        rg = max(abs(imag(str)),0.01*TC) * 2. 
        trr_values = range(real(str) - rg, stop = real(str) + rg, length=51) # we want the length to be an odd number
        tri_values = range(imag(str) - rg, stop = imag(str) + rg, length=51)
    else
        # println("I'm in regime (3), because ", abs(imag(saddle.tr)) )
        N = 51
        trr_values = range(real(tr_cd.min), stop=real(tr_cd.max), length= N )
        tri_values = range(min(imag(tr_cd.min),0.),stop = imag(tr_cd.max),length= N)  
    end
    
    relevant = false

    ti = saddle.ti
    
    S_values = [S(b, Ip, ti, tr, saddle.q) for tr in (trr_values' .+ im*tri_values)]   
    S_saddle = (S(b, Ip, saddle))

    # original integration domain (line)
    intdomain_tr = Curve2([
            (real(saddle.ti)/TC, -0.0),
            ( max( real(tr_cd.max)./TC, trr_values[end]/TC ) , -0.0)])
    
    con_S_saddle_real = Contour.contour(trr_values./TC, tri_values./TC, real.(S_values'), real.(S_saddle))
    for curve in lines(con_S_saddle_real)
            
        ip = intersection(curve, intdomain_tr)
        # println("all ips: $ip")
        
        for ipx in ip
            tr_ip = (ipx.x + ipx.y*im)*TC
            S_ip = S(b, Ip, ti, tr_ip, saddle.q)
            
            ascending = round(imag(S_ip - S_saddle),digits=5) >= 0
            crossing = crosses_point(curve,Point(reim(saddle.tr)./TC...), max(tri_values[2]-tri_values[1], 1.))
            # println(" tr check: asc $(ascending), crossing $crossing")

            if ascending && crossing
                relevant = true
                return relevant
            end
        end        

    end

    return relevant

end;


function plot_action_contours(beam::Beam,   
    s::Saddle,    
    ti_cd::ComplexDomain, tr_cd::ComplexDomain,
    ; all_saddles::Vector{Saddle}=Vector{Saddle}()
    )
    
    TC = TCycle(beam)
    q = s.q 
    Ntimes = 75

    println("consider the saddle point: ", s/TC)
    
    ### ionisation time
    plt_ti = let 
        tir_values = range(real(ti_cd.min)-0.1*TC, stop=real(ti_cd.max)+0.1*TC, length=Ntimes)
        tii_values = range(min(imag(ti_cd.min), 0.), stop = imag(ti_cd.max), length=Ntimes)   
        tr = s.tr

        Sv_values = [S_v(beam, Ip, ti, tr) for ti in (tir_values' .+ im*tii_values)]
        SvDrv_values = [dSv_dti(beam, Ip, ti, tr) for ti in (tir_values' .+ im*tii_values ) ]

        ImSv_values, _ = normalize(imag.(Sv_values))
        ReSv_values, _ = normalize(real.(Sv_values))

        plt1 = Plots.contour(tir_values./TC, tii_values./TC, ImSv_values
            , fill = true, aspectratio = 1, legend = false
            , levels = 30, linestyle = :solid
            , xlabel = "real ionisation times / Tcycle" , ylabel = "imag ionisation times / Tcycle"
            , title = "imag(Sv)");

        Plots.contour!(plt1, tir_values./TC, tii_values./TC, ReSv_values
            , fill = false, linestyle = :solid, c= :black, levels = 51)   

        Sv_saddle = (S_v(beam, Ip, s))

        mycurve = Curve2([(minimum(tir_values)./TC, 0.), (maximum(tir_values)./TC, 0.)])    

        con_S_saddle_real = Contour.contour(tir_values./TC, tii_values./TC, real.(Sv_values'), real.(Sv_saddle))

        for curve in lines(con_S_saddle_real)
            ip = intersection(curve,mycurve)
            for ipx in ip
                ti = (ipx.x + ipx.y*im)*TC
                Sv_ip = S_v(beam, Ip, ti, tr)
                ascending = imag(Sv_ip - Sv_saddle) > 0 
                
                crossthresh = max(tii_values[2]-tii_values[1], 1.5)
                crossespoint = crosses_point(curve,Point(reim(s.ti)./TC...), crossthresh )
                println("ti check: asc $(ascending), crossing $crossespoint => ti ", ascending && crossespoint )
            end
            plot!(plt1, coordinates(curve),  lw = 3)
        end 
        plot!(plt1, coordinates(mycurve), c = :red, lw = 3)
        
        ### plot all other saddle points as well
        scatter!(plt1,[reim(s.ti)./TC for s in all_saddles])
    end;
    
    ### recombination time
    
    plt_tr = let
        relevant = false

        str = s.tr
        rg = max(abs(imag(s.tr)),0.01*TC ) * 2. #
        trr_values = range(real(str) - rg, stop = real(str) + rg, length=51)
        tri_values = range(imag(str) - rg, stop = imag(str) + rg, length=51)

        #         println("step tri: ", tri_values[2]-tri_values[1])
            ti = s.ti

            S_values = [S(beam, Ip, ti, tr, q) for tr in (trr_values' .+ im*tri_values)]

            ImS_values, maxi = normalize(imag.(S_values))
            ReS_values, maxiR = normalize(real.(S_values))

            plt2 = Plots.contour(trr_values./TC, tri_values./TC, ImS_values
                , fill = true, aspectratio = 1, legend = false, levels = 30, linestyle = :solid
                , xlabel = "real return times / Tcycle" , ylabel = "imag return times / Tcycle", title = "imag(S)")

            Plots.contour!(trr_values./TC, tri_values./TC, ReS_values
                , fill = false, linestyle = :solid, c= :black, levels = 51)   

            S_saddle = (S(beam, Ip, s))
            con_S_saddle = Contour.contour(trr_values./TC, tri_values./TC, ImS_values', imag.(S_saddle)/maxi)
        
            # original integration domain (line)
            intdomain_tr = Curve2([
                    (real(s.ti)/TC, -0.0),
                    ( max( real(tr_cd.max)./TC, trr_values[end]/TC ) , -0.0)])
            @show intdomain_tr
        #         println("S saddle: ", S_saddle)
        con_S_saddle_real = Contour.contour(trr_values./TC, tri_values./TC, real.(S_values'), real.(S_saddle))
        for curve in lines(con_S_saddle_real)
            ip = intersection(curve, intdomain_tr)
            # println("all ips: $ip")
            for ipx in ip
                tr = round((ipx.x + ipx.y*im)*TC, digits=3)
                S_ip = S(beam, Ip, ti, tr, s.q)
                
                # println("ip: ", ipx)
                # println(" S ip: ", S_ip)
                # println(" diff: ", round(imag(S_ip - S_saddle), digits=4))
                
                ascending = round(imag(S_ip - S_saddle),digits=4) >= 0
                crossthresh =  max(tri_values[2]-tri_values[1], 1.)
                crossespoint = crosses_point(curve, Point(reim(s.tr)./TC...), crossthresh)
                println(" tr check: asc $(ascending), crossing $crossespoint")
                if ascending && crossespoint
                    relevant = true 
                end 
            end
           
            plot!(coordinates(curve), lw = 3)
        end
        println("==> tr relevant: ", relevant)
        plot!(coordinates(intdomain_tr), c = :red, lw = 3)
    #            plot!(xlims=(-0.25,-0.15), ylims = (-0.05,0.05))
        
        scatter!(plt2,[reim(s.tr)./TC for s in all_saddles])
        plt2
    end
    
    return plot(plt_ti, plt_tr, layout=(1,2), size =(800,300))
end