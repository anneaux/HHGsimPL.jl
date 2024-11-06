import ColorSchemes.turbo
function harmonic_color(q::T1, orders::Vector{T2}) where {T1 <: Real, T2 <:Real}
    if (maximum(orders)-minimum(orders)) == 0.
        return :blue 
    else 
        index = (q .- minimum(orders)) ./ (maximum(orders)-minimum(orders))
        return turbo[index]
    end
end


function drawzerolines(color::Symbol=:grey80)
  hline!([0],lc = color, label = "", z_order = :back)
  vline!([0],lc = color, label = "", z_order = :back)
end;


#### saddle points in the complex plane
function plot_sp_times_cp(saddlesTup::Vector{Tuple{Saddle,Bool}},
    TC::Float64=1.,
    thcs::Vector{Thc}=Vector{Thc}(),
    ti_cd::ComplexDomain = ComplexDomain(),
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(lambda=lambda)),
    tt_cd::ComplexDomain = ComplexDomain())
    
    plt_tt = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "travel time")
    plt_ti = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "ionization time")
    plt_tr = plot(legend = false, xlabel = "real part", ylabel = "imaginary part", title = "recombination time")
        
    for tup in saddlesTup
        s = tup[1]
        relevant = tup[2]
        
        col = harmonic_color(real(s.q), [real(tup[1].q) for tup in saddlesTup])
        # col = harmonic_color(Float64(real(s.q)), collect(15:40.))
        
        spmarker = (:o,2, relevant ? col : :white,stroke(2,col))
        
        ttr = real(s.tr-s.ti)./TC
        tti = imag(s.tr-s.ti)./TC
        
        ### travel time
        scatter!(plt_tt,(ttr,tti), aspectratio = 1., marker = spmarker)
        scatter!(plt_ti,reim(s.ti)./TC, aspectratio = 1., marker = spmarker)
        scatter!(plt_tr,reim(s.tr)./TC, aspectratio = 1., marker = spmarker)
    end

    plot_shape(plt, cd::ComplexDomain) =
        plot!(plt, Shape(real.([cd.min,cd.max,cd.max,cd.min])./TC,
            imag.([cd.min,cd.min,cd.max,cd.max])./TC), 
        fillalpha = 0,
        fillcolor = :yellow, 
        linecolor = :grey,
        linewidth = 3 )
    
    
    if !iszero(ti_cd.min)
        plot_shape(plt_ti, ti_cd)
        plot_shape(plt_tr, tr_cd)
    end
    
    thcmarker = (:diamond,3,:black,stroke(3,:black))
    for thc in thcs
        ttr = real(thc.trhc-thc.tihc)./TC
        tti = imag(thc.trhc-thc.tihc)./TC

        scatter!(plt_tt,(ttr,tti), aspectratio = 1., marker = thcmarker)
        scatter!(plt_ti,reim(thc.tihc)./TC, aspectratio = 1., marker = thcmarker)
        scatter!(plt_tr,reim(thc.trhc)./TC, aspectratio = 1., marker = thcmarker)
    end 
    
    return plot(plt_ti, plt_tr, plt_tt, layout=(1,3), size = (900,250))
end;


function plot_sp_times_cp(b::Beam, saddles::Vector{Saddle},
    thcs::Vector{Thc}=Vector{Thc}(),
    ti_cd::ComplexDomain = ComplexDomain(),
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
    tt_cd::ComplexDomain = ComplexDomain())

    TC = TCyle(b)
    saddleTups = [(s,true) for s in saddles]

    return plot_sp_times_cp(saddleTups, TC, thcs, ti_cd, tr_cd, tt_cd)
end


function plot_sp_times_cp(b::Beam, saddles_dict::Dict{T,Vector{Saddle}},
        thcs::Vector{Thc}=Vector{Thc}(),
        ti_cd::ComplexDomain = ComplexDomain(),
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        tt_cd::ComplexDomain = ComplexDomain()) where T <: Number

    plot_sp_times_cp(b, vcat(values(saddles_dict)...), thcs, ti_cd, tr_cd, tt_cd)
end

function plot_sp_times_cp(b::Beam, saddlesTup::Vector{Tuple{Saddle,Bool}},
    thcs::Vector{Thc},
    ti_cd::ComplexDomain = ComplexDomain(),
    tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
    tt_cd::ComplexDomain = ComplexDomain()) 
    
    return plot_sp_times_cp(saddlesTup, TCycle(beam), thcs, ti_cd, tr_cd, tt_cd)
end;


function plot_sp_times_cp(b::Beam, saddlesTup_dict::Dict{T,Vector{Tuple{Saddle,Bool}}},
        thcs::Vector{Thc}=Vector{Thc}(),
        ti_cd::ComplexDomain = ComplexDomain(),
        tr_cd::ComplexDomain = ComplexDomain(real(ti_cd.min)-imag(ti_cd.max)*im,ti_cd.max + TCycle(b)),
        tt_cd::ComplexDomain = ComplexDomain()) where T <: Number

    plot_sp_times_cp(b, vcat(values(saddlesTup_dict)...), thcs, ti_cd, tr_cd, tt_cd)
end


######## plot saddle points over real time (bit Anka style)
function plot_sp_times_real(saddlesTuples::Vector{Tuple{Saddle,Bool}}, TC::Float64=1.)   
    orders = [real(sTup[1].q) for sTup in saddlesTuples]
    
    ti_min = minimum(sTup -> real(sTup[1].ti), saddlesTuples)/TC
    ti_max = maximum(sTup -> real(sTup[1].ti), saddlesTuples)/TC
    tr_min = minimum(sTup -> real(sTup[1].tr), saddlesTuples)/TC
    tr_max = maximum(sTup -> real(sTup[1].tr), saddlesTuples)/TC 
    
    ioni_lims = (ti_min-0.1, ti_max+0.1)
    reco_lims = (tr_min-0.1, tr_max+0.1)

    plt_ioni = plot(legend = false
        , xlabel = "ionisation time / TC", ylabel = "harmonic order", xlims = ioni_lims)
    plt_reco = plot(legend = false
        , xlabel = "recombination time / TC", ylabel = "harmonic order", xlims = reco_lims)

    for tup in saddlesTuples
        s = tup[1]
        col = harmonic_color(real(s.q), orders)
        spmarker_i = (:diamond, 3,  tup[2] ? col : :white, tup[2] ? 1. : 1., stroke(1,col, 1.))
        spmarker_r = (:o, 3, tup[2] ? col : :white, tup[2] ? 1. : 1., stroke(1,col, 1.))
        
        scatter!(plt_ioni,[(real(s.ti)./TC, real(s.q))],
            marker = spmarker_i, label = false)
        scatter!(plt_reco,(real(s.tr)./TC, real(s.q))
            , marker = spmarker_r, label = false)
    end

    return plot(plt_ioni, plt_reco, layout = grid(1, 2, widths = [0.4,0.6 ]), size = (900, 400))

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

function plot_necklace_2D(necklace::Vector{LineSeg}, ti::ComplexF64, tr::ComplexF64)
    plt_imag = plot(xlabel = "Im(ti)", ylabel = "Im(tr)")
        plot!(plt_imag, [imag.((necklace[i].s.x, necklace[i].s.y)) for i in 1:length(necklace)],
            label = false, marker = true)
        scatter!(plt_imag, [imag.((ti,tr))], marker = (:red, :x, 5), label = "saddle point", legend = :bottom)

    plt_real = plot(xlabel = "Re(ti)", ylabel = "Re(tr)")
        plot!(plt_real, [real.((necklace[i].s.x, necklace[i].s.y)) for i in 1:length(necklace)],
            label = false, marker = true)
        scatter!(plt_real, [real.((ti,tr))], marker = (:red, :x, 5), label = false)

    return plot(plt_imag, plt_real, layout = (1,2), plot_title = "ti: $ti, tr: $tr", legend = :bottom)
end

