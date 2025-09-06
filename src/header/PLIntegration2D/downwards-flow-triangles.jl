#### OTHER UTILS
Base.getindex(z::Iterators.Zip, i) = (it -> getindex(it, i)).(z.is)
toR4(z1::ComplexF64, z2::ComplexF64) = [real(z1), imag(z1), real(z2), imag(z2)]
toR4(z1::Number, z2::Number) = [real(z1), imag(z1), real(z2), imag(z2)]




#### POINTS
mutable struct PointA{T}
    x::T
    y::T
    active::Bool
end
PointA(x,y) = PointA(x,y,true);

dist(p1::PointA, p2::PointA) = norm([p1.x-p2.x, p1.y-p2.y])

import Base.isequal
isequal(p1::PointA,p2::PointA) = isequal(p1.x,p2.x) && isequal(p1.y,p2.y)
Base.hash(p::PointA, h::UInt) = hash([p.x,p.y], h)

function float2complex(p::PointA)
    return PointA(complex(p.x), complex(p.y), p.active)
end
toR4(p::PointA) = toR4(p.x, p.y)




#### TRIANGLES
mutable struct TriangleA
    coord::MVector{3,Int} ### this could be an MVector
    active::Bool
    TriangleA(coord::AbstractVector{T}) where T <:Integer = new(coord, true)
end

function project_onto_triangle(base::AbstractVector, points::AbstractVector) # what are these types?
	    @assert length(base) == 3 "Need exactly 3 vertices"

	    # Convert to R^4
	    vs = [toR4(z.x, z.y) for z in base]

	    # Basepoint
	    v1, v2, v3 = vs
	    e1 = v2 - v1
	    e2 = v3 - v1

	    # Gram-Schmidt to find orthonormal basis
	    u1 = e1 / norm(e1)
	    e2_proj = e2 - (u1 ⋅ e2) * u1
	    u2 = e2_proj / norm(e2_proj)

	    # Projection function
	    proj(v) = [(u1 ⋅ (v - v1)), (u2 ⋅ (v - v1))]

	    # Project all points
	    ps = [proj(toR4(p.x, p.y)) for p in points]

	    return ps
end

using MiniQhull
flags_original = "qhull d Qt Qbb Qc Qz"

function subdivide_triangle_new!(points, t_vertices::TriangleA; Δ::Real=0.5, delaunay_flags=flags_original)
    
	    vertices_here = Vector{Int}()
	    push!(vertices_here, t_vertices.coord...)
	    v1,v2,v3 = t_vertices.coord
	    t_edges = [(v1,v2),(v2,v3),(v3,v1)] ### indices of the specific one I'm looking at
	    for k in t_edges
	        p1,p2 = points[[k...]]
	        d = dist(p1,p2)
	        if d > Δ ### subdivide
	            
	            ### how many points to insert in between?
	            ### wanna make sure that the new distance is not larger than Δ
	            n = floor(Int64, d/Δ)
	            xvals = range(p1.x, stop = p2.x, length = n+2)
	            yvals = range(p1.y, stop = p2.y, length = n+2)
	            new_points = [PointA(c...) for c in zip(xvals[2:end-1], yvals[2:end-1])] # I don't want to include the actual points themselves

	            for np in new_points
	                already_created = findall(p->isequal(p, np), points) 
	                if isempty(already_created)
	                    push!(points, np)
	                    push!(vertices_here, length(points))
	#                     println("point is new. new index: $(length(points)),  length points: $(length(points))")
	                else
	#                         println("point is already created. length points: $(length(points))")
	                    np_vertex = already_created[1]
	                    push!(vertices_here, np_vertex)
	                    ### add a pointer to the existing point
	                end
	            end

	        end
	    end
	    
	    if length(vertices_here) == 3
	        return points, Vector{TriangleA}[]
	    else
	        ### project only the points that are relevant for this new triangle
	        projected_points = project_onto_triangle(points[t_vertices.coord], points[vertices_here])
	        connections_here = delaunay(projected_points, flags_original) ### they are in the index system of this triangle!

	        ### check that they are not degenerate!	        
	        connections_here_vectors = [ collect(col) for col in eachcol(connections_here)]
	        filter!(ind -> triangle_area_new(projected_points[ind]) > 1e-5,
	            connections_here_vectors)

	        new_triangles = [TriangleA(vertices_here[inds]) for inds in connections_here_vectors]
	        return points, new_triangles
	    end
end


function dissect_thimbles(triangles)
    active_triangles = filter(s->s.active, triangles)
    thimbles = Vector()
    visited = falses(length(active_triangles))
   
    function dfs!(active_triangles, visited, i_start)
        stack = [i_start] # "open ends" to explore (= the four corner point indices of the quadrilateral)
        trace = Int[] # trace is more like a carpet now I guess
        while !isempty(stack)
            v = pop!(stack)

            if !visited[v]
			# #                 @show v
                visited[v] = true
                push!(trace, v)

                ### find the neighbouring quads to the vertices
                for vertex in 1:3
                    # find all quads that corner to you
                    nexts = findall(sim -> in(active_triangles[v].coord[vertex], sim.coord), active_triangles)
					#                     @show nexts
                    append!(stack, filter(!isequal(v),nexts)) # obviously ignore the quad that you're in atm
                end
                unique!(stack)
            end
        end
        return trace
    end
    
    for i in 1:length(active_triangles)
		#     i = 1
        if !visited[i]
            trace = dfs!(active_triangles, visited, i)
            push!(thimbles, active_triangles[trace])
        end
    end
    
    return thimbles
end


function triangle_area_new(points::AbstractVector)
    @assert length(points) == 3
    abs((points[2][1]-points[1][1])*(points[3][2]-points[1][2]) - (points[2][2]-points[1][2])*(points[3][1]-points[1][1])) / 2
end
# triangle_area(points::AbstractVector{PointA}) = triangle_area([xy(p) for p in points])


function integrate_triangle(points, triangle, integrand; order=1,dim=2)
    int = zeros(ComplexF64,dim)
    area = triangle_area_new(xy.(points[triangle.coord]))
    wi = area/3

    v1,v2,v3 = triangle.coord
    t_edges = [(v1,v2),(v2,v3),(v3,v1)] ### indices of the specific one I'm looking at
    for k in t_edges
        p1,p2 = points[[k...]]
        ### is this fine or does this need to be done in the projection???
        midpoint = ((p1.x+p2.x)/2, (p1.y+p2.y)/2) 
        int .+= wi*integrand(midpoint...)
    end    
    return int
end








#### DOWNWARDS FLOW


function gradN(
    f_grad::Function,
    ti::ComplexF64, tr::ComplexF64,
    thresh::Float64 = 1.)

    g = f_grad(ti,tr)
    if norm(g) > thresh # bit lower than the gradient at the saddle point
        return LinearAlgebra.normalize(g)
    else 
        return g
    end
end;



function initialise_time_grid(timin, timax, ttmin, ttmax, Δ;
        flow_bounds=[true, true, true, true])
    
    points = [
        PointA(timin, timin+ttmin, flow_bounds[1]), 
        PointA(timin, timin+ttmax, flow_bounds[2]), 
        PointA(timax, timax+ttmax, flow_bounds[3]),
        PointA(timax, timax+ttmin, flow_bounds[4])]      
    
    pts = [SVector(p.x, p.y) for p in points]
    connections = delaunay(pts) #, "qhull d Qbb Qc QJ Pp")
    #     "qhull d Qbb Qc QJ Pp"
	#         display(connections)
    triangles = [TriangleA(inds) for inds in eachcol(connections)]
	#     @show points
    points = [float2complex(p) for p in points]
	#     @show points
    
    ### filter out degenerate points!
    n_old = length(triangles)
    n_new = n_old + 1

    i_while = 0
    while (n_old != n_new)
        n_old = n_new 
        i_while += 1
        i = 1    
    
        for i_t in eachindex(triangles)
            triangle = triangles[i_t]
            if triangle.active
		#                 points_this_triangle = points[triangle.coord]
                points, new_connections = subdivide_triangle_new!(points, triangle, Δ=Δ)
                if !isempty(new_connections)
                    triangle.active = false
                    append!(triangles, new_connections)
                end

            end
            if i>30 break end
        end

        if i_while > 10 break end
        filter!(sim->sim.active, triangles)
        n_new = length(triangles)

    end         
    return points, triangles
end



function initialise_triangles_on_real_plane(xmin, xmax, ymin, ymax, Δx, Δy)
    xrange = range(real(xmin), stop = real(xmax), length = Int(ceil(real(xmax-xmin)/Δx)))
    yrange = range(real(ymin), stop = real(ymax), length = Int(ceil(real(ymax-ymin)/Δy)))

    points_ini = Vector{PointA}()
    for xr in collect(xrange)
        for yr in collect(yrange)
            push!(points_ini, PointA(xr, yr,true))
        end
    end
    pts = [SVector(p.x, p.y) for p in points_ini]
    connections = delaunay(pts) #, "qhull d Qbb Qc QJ Pp")
    #     "qhull d Qbb Qc QJ Pp"
	#         display(connections)

	### filter out degenerate points!



    return points_ini, connections
end


function flow_down!(triangles, points::Vector{PointA{T}},
        f::Function,
        f_grad::Function;
        threshold::Real=0.5, # for normalisation of thr gradient
        δ::Real=0.5, # flowstepfactor
        h_threshold::Real=-20.
        ) where T<:Number

    for i1 in 1:length(points)
        if points[i1].active # for the active points
            step = -δ .* gradN((ti,tr) -> conj.(complex.(f_grad(ti,tr))), points[i1].x +0im, points[i1].y +0im, threshold)
            points[i1].x += step[1]
            points[i1].y += step[2]
        end
    end

    for i2 in eachindex(triangles)
        if triangles[i2].active
            for v in triangles[i2].coord
                if real(f(points[v].x, points[v].y)) < h_threshold 
                    triangles[i2].active = false # am I sure that I want to turn the whole simplex inactive?
                    points[v].active = false
                end
            end
        end
    end
end


function subdivide_triangles!(points, triangles, subdividethreshold)
    n_old = length(triangles)
    n_new = n_old + 1

    i_while = 0
    while (n_old != n_new)
        n_old = n_new 
        i_while += 1
        i = 1

        ### in subdivide sub routine ideally
        for i_t in eachindex(triangles)
            triangle = triangles[i_t]
            if triangle.active
                if all([p.active for p in points[triangle.coord]])
                    points, new_connections = subdivide_triangle_new!(points, triangle, Δ=subdividethreshold)
					# n_actives = sum([t.active for t in triangles])
                    if !isempty(new_connections)
                        triangle.active = false
                        append!(triangles, new_connections)
                    end
                elseif !any([p.active for p in points[triangle.coord]])
                    triangle.active = false
                    println("here i have turnt a triangle inactive")
                end
            end
        end
        
		#  if i_while > 10 break end
        filter!(sim->sim.active, triangles)
        n_new = length(triangles)
    end
end