import ColorSchemes.turbo
function harmonic_color(q::T, orders::Vector{T}) where T <: Real
    if (maximum(orders)-minimum(orders)) == 0.
        return :blue 
    else 
        index = (q .- minimum(orders)) ./ (maximum(orders)-minimum(orders))
        return turbo[index]
    end
end




function plot_sp_times_cp(b::Beam, saddles::Vector{Saddle}, thcs::Vector{Thc}=Vector{Thc}(),
        ti_cd::ComplexDomain = ComplexDomain(),
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        tt_cd::ComplexDomain = ComplexDomain()) 
    
    plt_tt = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "travel time")
    plt_ti = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "ionization time")
    plt_tr = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "recombination time")
    
    for s in saddles
        col = harmonic_color(real(s.q), [real(s.q) for s in saddles])
        spmarker = (:o,2,col,stroke(2,col))
        
        ttr = real(s.tr-s.ti)./TCycle(b)
        tti = imag(s.tr-s.ti)./TCycle(b)
        
        ### travel time
        scatter!(plt_tt,(ttr,tti), aspectratio = 1., marker = spmarker)
        scatter!(plt_ti,reim(s.ti)./TCycle(b), aspectratio = 1., marker = spmarker)
        scatter!(plt_tr,reim(s.tr)./TCycle(b), aspectratio = 1., marker = spmarker)
    end

    plot_shape(plt, cd::ComplexDomain) =
        plot!(plt, Shape(real.([cd.min,cd.max,cd.max,cd.min])./TCycle(b),
            imag.([cd.min,cd.min,cd.max,cd.max])./TCycle(b)), 
        fillalpha = 0,
        fillcolor = :yellow, 
        linecolor = :grey,
        linewidth = 3 )
    
    plot_shape(plt_ti,ti_cd)
    
    if !iszero(ti_cd.min)
        plot_shape(plt_tr, tr_cd)
    end
    
    thcmarker = (:diamond,3,:black,stroke(3,:black))
    for thc in thcs
        ttr = real(thc.trhc-thc.tihc)./TCycle(b)
        tti = imag(thc.trhc-thc.tihc)./TCycle(b)

        scatter!(plt_tt,(ttr,tti), aspectratio = 1., marker = thcmarker)
        scatter!(plt_ti,reim(thc.tihc)./TCycle(b), aspectratio = 1., marker = thcmarker)
        scatter!(plt_tr,reim(thc.trhc)./TCycle(b), aspectratio = 1., marker = thcmarker)
    end 
    
    return plot(plt_ti, plt_tr, plt_tt, layout=(1,3), size = (900,250))
end;

function plot_sp_times_cp(b::Beam, saddles_dict::Dict{T,Vector{Saddle}},
        thcs::Vector{Thc}=Vector{Thc}(),
        ti_cd::ComplexDomain = ComplexDomain(),
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        tt_cd::ComplexDomain = ComplexDomain()) where T <: Number

    plot_sp_times_cp(b, vcat(values(saddles_dict)...), thcs, ti_cd, tr_cd, tt_cd)
end




function plot_necklace_3D(necklace::Vector{LineSeg}, 
    # b::Beam, Ip::Float64,
    q::Number,
    # ti::ComplexF64, tr::ComplexF64,
    ti_cd::ComplexDomain, tr_cd::ComplexDomain,
    contourline::Curve2{Tuple{T, T}} where T<:Real
    ; Ntimes::Int64=100)
    
    tir_values = realrange(ti_cd, Ntimes); tii_values = imagrange(ti_cd, Ntimes)
    trr_values = realrange(tr_cd, Ntimes); tri_values = imagrange(tr_cd, Ntimes)    

    plt_1 = plot(xlabel = "Re(ti)", ylabel = "Re(tr)", zlabel = "Im(ti)")
    plt_2 = plot(xlabel = "Re(ti)", ylabel = "Re(tr)", zlabel = "Im(tr)")
    
    ### original integration (half) plane
    integrationplane(ti, tr) = (real(tr) > real(ti)) ? 0. : missing
    surface!(plt_1, tir_values, trr_values, 
        integrationplane.(tir_values', trr_values), label = "integration plane", 
        alpha = 0.5, color = :red, cbar = false)
    surface!(plt_2, tir_values, trr_values, 
        integrationplane.(tir_values', trr_values), label = "integration plane", 
        alpha = 0.5, color = :red, cbar = false)
    
    ### level line for the saddle point
    xs, ys = coordinates(contourline)
    plot!(plt_1, xs, ys, zeros(length(xs)), color = :blue, linewidth = 2, label = "Im(S) contour line at the SP")
    plot!(plt_2, xs, ys, zeros(length(xs)), color = :blue, linewidth = 2, label = "Im(S) contour line at the SP")
    ### necklace
    qcol = harmonic_color(Float64(q), collect(15.:40))
    mymarker = (:o,1,qcol,stroke(0.5,:black))
    scatter3d!(plt_1,[(real(necklace[n].s.x), real(necklace[n].s.y), imag(necklace[n].s.x)) 
            for n in 1:length(necklace)], marker = mymarker, label = "loop")
    scatter3d!(plt_2,[(real(necklace[n].s.x), real(necklace[n].s.y), imag(necklace[n].s.y)) 
            for n in 1:length(necklace)], marker = mymarker, label = "loop")
   
    plot!(plt_1, zlims = (-2,:auto), title = "view 1")    
    plot!(plt_2, title = "view 2")
    plot!(plt_1, xguide = "Re(ti)", ylabel = "Re(tr)", zlabel = "Im(ti)")

    display(plt_1); display(plt_2)
#     plot(plt_1, plt_2, layout = (1,2), size = (800,400)) 
    
end