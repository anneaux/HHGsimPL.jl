### for GLMakie
### getting points that I can plot "transparent" just to create the correct axis dimensions
function get_init_points(xlims::Tuple{Real,Real}, ylims::Tuple{Real,Real})
    [
        (xlims[1],ylims[1]),
        (xlims[2],ylims[2])
    ]
end

function get_init_points(xlims::Tuple{Real,Real},
    ylims::Tuple{Real,Real},
    zlims::Tuple{Real,Real},
    scaling=(1., 1., 1.))

    [
        (xlims[1],ylims[1],zlims[1]) ./scaling,
        (xlims[2],ylims[2],zlims[2]) ./scaling
    ]

end

using Printf

function get_ticks(lims::Vector, scaling::Union{Tuple,Vector}, N::Vector{Int64}=[5,5,5] )
    xticks = range(lims[1][1], stop = lims[1][2], length = N[1]) 
    yticks = range(lims[2][1], stop = lims[2][2], length = N[2]) 
    zticks = range(lims[3][1], stop = lims[3][2], length = N[3]) 
    #     @show yticks
    ticks = (xticks, yticks, zticks ) ./ scaling

    xlabs = [@sprintf("%.2f",x) for x in xticks]
    ylabs = [@sprintf("%.2f",y) for y in yticks]
    zlabs = [@sprintf("%.2f",z) for z in zticks]
    
    labels = (xlabs, ylabs, zlabs )
    
    return (ticks, labels)
end



#### 
function get_ticks(lims::Vector;
        xreverse::Bool=false, yreverse::Bool=false,zreverse::Bool=false,
        scaling::Union{Tuple,Vector}=ones(3),
        Nticks::Vector{Int64}=[5,5,5]
    )
    xticks = collect(range(lims[1][1], stop = lims[1][2], length = Nticks[1]) )
    yticks = collect(range(lims[2][1], stop = lims[2][2], length = Nticks[2]) )
    zticks = collect(range(lims[3][1], stop = lims[3][2], length = Nticks[3]) )
    ticks = ( xreverse ? lims[1][1] .- xticks : xticks,
        yreverse ? lims[2][1] .- yticks : yticks,
        zreverse ? lims[3][1] .- zticks : zticks) ./ scaling

    xlabs = [@sprintf("%.2f",x) for x in xticks]
    ylabs = [@sprintf("%.2f",y) for y in yticks]
    zlabs = [@sprintf("%.2f",z) for z in zticks]
    
    labels = (xlabs, ylabs, zlabs )
    
    return (Tuple(ticks), labels) 
end

function get_init_points(xlims::Tuple{Real,Real},
    ylims::Tuple{Real,Real},
    zlims::Tuple{Real,Real};
    xreverse::Bool=false,yreverse::Bool=false,zreverse::Bool=false,
    scaling::Union{Tuple,Vector}=ones(3))
    
    return [
        (xreverse ? 0. : xlims[1],
        yreverse ? 0. : ylims[1],
        zreverse ? 0. : zlims[1]) ./scaling,
        (xreverse ? xlims[1] - xlims[2] : xlims[2],
        yreverse ? ylims[1] - ylims[2] : ylims[2],
        zreverse ? zlims[1] - zlims[2] : zlims[2] ) ./scaling
    ]  
end

function initialise_axes!(ax::LScene, xlims::Tuple, ylims::Tuple, zlims::Tuple;
        xlabel::String="x", ylabel::String="y", zlabel::String="z",
        kwargs...
    #         xreverse::Bool=true,
    #         scaling::Tuple=(1.,1., 1.),
    #         Nticks::Vector{Int64}=[5,5,5]
    )
    
    axis = ax.scene[OldAxis]
    axis[:names, :axisnames] = (xlabel, ylabel, zlabel)

    ticks = get_ticks([xlims, ylims, zlims]; kwargs...)
    ### for when I'm on my office Pc
    # axis[:ticks][][:ranges] = ticks[1]
    # axis[:ticks][][:labels] = ticks[2]

    ### for when I'm on my laptop
    axis[:ticks][:ranges] = ticks[1]
    axis[:ticks][:labels] = ticks[2]
   
    ### adding some initial point just to get the axis limits right...
    points_init = get_init_points(xlims, ylims, zlims;  kwargs...)
    scatter!(ax, [Point3f(p) for p in points_init], color = :transparent)
end


### my typical coordinate projections for 2d and 3d plots

coords_2d = [
    p -> (real(p[1]), imag(p[1])),
    p -> (real(p[2]), imag(p[2])),  
    p -> (real(p[1]), real(p[2])),  
    p -> (imag(p[1]), imag(p[2])) 
];


coords_3d = [
    p -> (real(p[1]), real(p[2]), imag(p[1])),
    p -> (real(p[2]), real(p[1]), imag(p[2])),
];

nothing
